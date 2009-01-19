#include "hermes2d.h"
#include <GL/glut.h>
#include <unistd.h>
#include <getopt.h>


template<class Base>
class LinView : public Base
{
public:
  
  LinView(int width, int height, int n, char** names) 
    : current(0), num(n), names(names),
      Base("LinView", 0, 0, width, height) {}
  
  void switch_to(int n)
  {
    Base::load_data(names[n]);
    update_title();
  }
  
  int get_num() const { return num; }
    
protected:
  
  int num, current;
  char** names;

  virtual void on_special_key(int key, int x, int y)
  {
    switch (key)
    {
      case GLUT_KEY_PAGE_UP:
        if (current > 0) switch_to(--current);
        break;
        
      case GLUT_KEY_PAGE_DOWN:
        if (current < num-1) switch_to(++current);
        break;

      case GLUT_KEY_HOME:
        switch_to(current = 0);
        break;

      case GLUT_KEY_END:
        switch_to(current = num-1);
        break;
        
      default: Base::on_special_key(key, x, y);
    }
  }
  
  void update_title()
  {
    char buffer[1000];
    sprintf(buffer, "LinView - %s", names[current]);
    Base::set_title_internal(buffer);
  }

};


void print_help()
{
  printf("Usage: linview [OPTIONS] FILE [FILE...]\n"
         "Displays saved Hermes2D Linearizer files.\n\n"
         "  -w, --width=pixels   initial window width\n"
         "  -h, --height=pixels  initial window height\n"
         "  -s, --scale          display scale\n"
         "  -m, --min=value      scale minimum value\n"
         "  -M, --max=value      scale maximum value\n"
         "  -v, --video          capture video frames on keypress\n"
         "  -f, --filename       video frame file name template\n"
         "  -q, --highquality    render high quality video\n"
         "  -t, --scale-format   scale printf number format (default %%.3g)\n"
         "  -p, --pal-steps=val  number of palette steps\n"
         #ifdef CONTOURS
         "  -k, --cont-steps=val  interval of contour steps\n"
         "  -o, --cont-orig=val   origin of contours\n"
         #endif
         "\n" );

}


bool exists(const char* filename)
{
  FILE* f = fopen(filename, "r");
  if (f == NULL) return false;
  fclose(f);
  return true;
}


int main(int argc, char* argv[])
{
  verbose_mode = false;

  struct 
  {
    int width, height, palsteps;
    double consteps, conorig;
    bool scale, video, hq, contours;
    double min, max;
    char filename[256];
    char scaleformat[50];
  }
  options;
  
  options.width = 1000;
  options.height = 750;
  options.scale = false;
  options.video = false;
  options.min = 0.0;
  options.max = 0.0;
  options.hq = false;
  options.palsteps = 50;
  options.contours = false;
  options.consteps = 0.1;
  options.conorig = 0.0;
  strcpy(options.filename, "frame%04d.bmp");
  strcpy(options.scaleformat, "%.3g");
  
  static const option long_opts[] =
  {
    { "width",        1, NULL, 'w' },
    { "height",       1, NULL, 'h' },
    { "min",          1, NULL, 'm' },
    { "max",          1, NULL, 'M' },
    { "scale",        0, NULL, 's' },
    { "video",        0, NULL, 'v' },
    { "help",         0, NULL, '?' },
    { "filename",     1, NULL, 'f' },
    { "highquality",  0, NULL, 'q' },
    { "scale-format", 1, NULL, 't' },
    { "pal-steps",    1, NULL, 'p' },
    { "cont-steps",    1, NULL, 'k' },
    { "cont-orig",     1, NULL, 'o' },
    { NULL,           0, NULL,  0  }
  };
  
  int opt, long_idx;
  while ((opt = getopt_long(argc, argv, "w:h:m:M:svqf:t:p:k:o:?", long_opts, &long_idx)) != -1)
  {
    switch (opt)
    {
      case 'w': options.width = atoi(optarg); break;
      case 'h': options.height = atoi(optarg); break;
      case 'm': sscanf(optarg, "%lf", &options.min); break;
      case 'M': sscanf(optarg, "%lf", &options.max); break;
      case 's': options.scale = true; break;
      case 'v': options.video = true; break;
      case 'f': strcpy(options.filename, optarg); break;
      case 'q': options.hq = true; break;
      case 't': strcpy(options.scaleformat, optarg); break;
      case 'p': options.palsteps = atoi(optarg); options.contours = true; break;
      case 'k': options.consteps = atof(optarg); options.contours = true; break;
      case 'o': options.conorig = atof(optarg); break;
      default:  print_help(); return EXIT_FAILURE;
    }
  }
  
  if (optind >= argc) 
  {
    printf("No files to process. To get help, run with --help\n");
    return EXIT_FAILURE;
  }
  
  LinView<WHATVIEW> linview(options.width, options.height, argc-optind, argv+optind);
  
  linview.show_scale(options.scale);
  linview.set_scale_format(options.scaleformat);
  if (options.min != options.max)
    linview.set_min_max_range(options.min, options.max);
  linview.set_num_palette_steps(options.palsteps);

  #ifdef CONTOURS
  if (options.contours)
    linview.show_contours(options.consteps, options.conorig);
  #endif

  linview.switch_to(0);

  if (options.video)
  {
    linview.wait_for_keypress();
    for (int i = 0; i < linview.get_num(); i++)
    {
      char filename[100];
      sprintf(filename, options.filename, i+1);
      if (!exists(filename))
      { 
        linview.switch_to(i);
        //usleep(100000); // bad intel GL driver
        linview.save_screenshot(filename, options.hq);
      }
      else
        info("%s already exists: skipping", filename);
    }
  }
  
  View::wait();
  return EXIT_SUCCESS;
}

