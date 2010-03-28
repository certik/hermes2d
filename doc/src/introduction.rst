============
Introduction
============

Hermes2D is a free C++/Python library for rapid prototyping of
adaptive FEM and *hp*-FEM solvers for partial differential equations (PDE),
developed by the `hp-FEM group <http://hpfem.org/>`_ at the University of 
Nevada, Reno. The library is available under the GPL license (Version 2, 1991).
Although Hermes is much younger than other FEM packages, it is loaded with 
unique technology and its user base is growing fast. In the following, 
we will abbreviate Hermes2D with Hermes. 

For those who speak other languages than C++, there is an **interactive 
GUI** `Agros2D <{http://hpfem.org/hermes2d/>`_. We also provide 
an **interactive online lab** (`accessible here <http://nb.femhub.org/>`_) where
you can compute with Hermes via any web browser without even installing it 
(the CPU time is on us). 

Prior to reading this document, we recommend that you install Hermes using instructions on 
its `home page <http://hpfem.org/hermes2d/>`_, and subscribe to the `mailing list 
<http://groups.google.com/group/hermes2d/>`_. Our mailing list is a very active place where 
you should get all answers quickly. 


The best way of reading this tutorial is to run the code at the same time. 
After making your way through the tutorial, you may want to view the directory 
`examples/ <http://hpfem.org/git/gitweb.cgi/hermes2d.git/tree/HEAD:/examples>`_ 
that contains a variety of different PDE models that may help you to get started with your own 
applications. If you create an interesting model using Hermes, let us know and we 
will be happy to add it to the existing examples. 

The source code can be 
viewed in the `git repository <http://hpfem.org/git/gitweb.cgi/hermes2d.git/tree>`_, 
and all tutorial examples can be found in the directory 
`tutorial/ <http://hpfem.org/git/gitweb.cgi/hermes2d.git/tree/HEAD:/tutorial>`_.
For the 1D and 3D codes, see the `Hermes1D <http://hpfem.org/hermes1d/>`_ and 
`Hermes3D <http://hpfem.org/hermes3d/>`_ home pages, respectively.

What is Hermes
--------------

Main strengths of Hermes are **adaptive hp-FEM methods**,
adaptivity for time-dependent problems on **dynamical hp-meshes**, 
monolithic discretization of arbitrary multiphysics problems via a novel **multimesh hp-FEM**,
and unprecedented **interactive web accessibility**. 
The following list gives more details so that you can decide whether Hermes 
may be the library that you are looking for: 

* **Mature hp-adaptivity algorithms**. Hermes puts a major emphasis on credibility of results, i.e., on error control and automatic adaptivity. Practitioners know well how painful it is to use automatic adaptivity in conjunction with standard lower-order approximations such as linear or quadratic elements - the error decreases somehow during a few initial adaptivity steps, but then it slows down and it does not help to invest more unknowns or CPU time. This is typical for low-order methods. In contrast to this, the exponentially-convergent adaptive *hp*-FEM and *hp*-DG do not have this problem - the error drops steadily and fast during adaptivity all the way to the desired accuracy. 

Typical convergence comparison of h-FEM with linear and quadratic elements, and the hp-FEM on a log-log scale:

.. image:: img/intro/conv_dof.png
   :align: center
   :width: 600
   :height: 400
   :alt: Typical convergence curves of FEM with linear and quadratic elements and hp-FEM

Same graphs as above but now in terms of CPU time:

.. image:: img/intro/conv_cpu.png
   :align: center
   :width: 600
   :height: 400
   :alt: CPU convergence graph.

* **Wide applicability**. Hermes is completely PDE-independent. Many FEM codes are designed to solve some narrow class of PDE problems (such as elliptic equations, fluid dynamics, electromagnetics etc.). In contrast to that, Hermes does not employ any technique or algorithm that would only work for some particular class of PDE problems. Automatic adaptivity is guided by a universal computational a-posteriori error estimate that works in the same way for any PDE. Of course this does not mean that it performs equally well on all PDE - some equations simply are more difficult to solve than others. However, Hermes allows you to tackle an arbitrary PDE or multiphysics PDE system. Visit the `hp-FEM group home page <http://hpfem.org/>`_ and especially the `gallery <http://hpfem.org/gallery/>`_ to see numerous examples.

.. image:: img/intro/ns.jpg
   :align: center
   :width: 650
   :height: 300
   :alt: Image of incompressible viscous flow.


* **Arbitrary-level hanging nodes**. Hermes has a unique original methodology for handling arbitrary-level hanging nodes. This means that extremely small elements can be adjacent to very large ones. When an element is refined, its neighbors are never split forcefully as in conventional adaptivity algorithms. It is well known that approximations with one-level hanging nodes are more efficient compared to regular meshes. However, the technique of arbitrary-level hanging nodes brings this to a perfection.

.. image:: img/intro/ord_2d_c.png
   :align: center
   :width: 370
   :height: 350
   :alt: Illustration of arbitrary-level hanging nodes.

.. ######
    .. image:: img/mixer-mesh.png
       :align: right
       :width: 300
       :height: 300
       :alt: Illustration of arbitrary-level hanging nodes.

    .. raw:: html

       <hr style="clear: both; visibility: hidden;">

* **Multimesh hp-FEM**. Various physical fields or solution components in multiphysics problems can be approximated on individual meshes, combining quality $H^1$, $H(curl)$, $H(div)$, and $L^2$ conforming higher-order elements. Due to a unique original methodology, no error is caused by operator splitting, transferring data between different meshes, and the like. The following figure illustrates a coupled problem of heat and moisture transfer in massive concrete walls of a nuclear reactor vessel. 

.. image:: img/intro/hm-sln-frame.png
   :align: left
   :width: 500
   :height: 410
   :alt: Illustration of multimesh hp-FEM.

.. image:: img/intro/hm-mesh-frame.png
   :align: right
   :width: 500
   :height: 410
   :alt: Illustration of multimesh hp-FEM.

.. raw:: html

   <hr style="clear: both; visibility: hidden;">

* **Dynamical meshes for time-dependent problems**. In time-dependent problems, different physical fields or solution components can be approximated on individual meshes that evolve in time independently of each other. Due to a unique original methodology, no error is caused by transfering solution data between different meshes and time levels. No such transfer takes place in the multimesh *hp*-FEM - the discretization of the time-dependent PDE system is monolithic. 

.. image:: img/intro/flame.jpg
   :align: center
   :width: 700
   :height: 360
   :alt: Adaptive hp-FEM with dynamical meshes for a flame propagation problem. 

* **Interactive web usage**. You can use Hermes remotely via any web browser, using our `interactive online lab <http://nb.femhub.org/>`_. Your hardware will not be used since the online lab is powered by the University of Nevada, Reno (UNR) high-performance computing facility (`Research Grid <http://hpc.unr.edu/wiki/index.php/Main_Page>`_). You can compute with Hermes using an iPhone if you like. Sound too good to be true? Try it. 

.. image:: img/intro/iphone_large.png
   :align: center
   :width: 250
   :height: 450
   :alt: Hermes in iPhone.

See the `Hermes home page <http://hpfem.org/main/hermes.php>`_ for more information. An overview of books, 
journal articles, conference proceedings papers and talks about Hermes and adaptive *hp*-FEM can be 
found in its `publications section <http://hpfem.org/publications/>`_.

Citing Hermes
-------------

If you use Hermes for your work, please be so kind to include some of the references below as appropriate.

Monograph:

::

    @Book{Hermes-book,
      author = {P. Solin, K. Segeth, I. Dolezel},
      title = {Higher-Order Finite Element Methods},
      publisher = {Chapman & Hall / CRC Press},
      year = {2003}
    }

Reference to the Hermes project:

::

    @Manual{Hermes-project,
      title =  {Hermes - Higher-Order Modular Finite Element System (User's Guide)},
      author = {P. Solin et al.},
      url =    {http://hpfem.org}
    }

Underlying algorithms (hanging nodes, adaptivity, shape functions):

::

    @Article{Hermes-hanging-nodes,
      author = {P. Solin, J. Cerveny, I. Dolezel},
      title = {Arbitrary-Level Hanging Nodes and Automatic Adaptivity in the hp-FEM},
      journal = {Math. Comput. Simul.},
      volume = {77},
      year = {2008},
      pages = {117 - 132},
      doi = {doi:10.1016/j.matcom.2007.02.011}
    }

::

    @Article{Hermes-adaptivity,
      author = {P. Solin, D. Andrs, J. Cerveny, M. Simko},
      title = {PDE-Independent Adaptive hp-FEM Based on Hierarchic Extension of Finite Element Spaces},
      journal = {J. Comput. Appl. Math.},
      status = {accepted},
      year = {2009},
    }

::

    @Article{Hermes-shape-functions,
      author = {P. Solin, T. Vejchodsky},
      title = {Higher-Order Finite Elements Based on Generalized Eigenfunctions of the Laplacian},
      journal = {Int. J. Numer. Methods Engrg},
      volume = {73},
      year = {2007},
      pages = {1374 - 1394}
    }

Topical papers from various application areas:

::

    @Article{Hermes-multiphysics,
      author = {P. Solin, L. Dubcova, J. Kruis},
      title = {Adaptive hp-FEM with Dynamical Meshes for Transient Heat and Moisture Transfer Problems},
      journal = {J. Comput. Appl. Math},
      doi = {doi 10.1016/j.cam.2009.07.025},
      year = {2009}
    }

::

    @Article{Hermes-solid-mechanics,
      author = {P. Solin, J. Cerveny, L. Dubcova, D. Andrs},
      title = {Monolithic Discretization of Linear Thermoelasticity Problems via Adaptive Multimesh hp-FEM},
      journal = {J. Comput. Appl. Math},
      doi = {doi 10.1016/j.cam.2009.08.092},
      year = {2009}
    }

::

    @Article{Hermes-electromagnetics,
      author = {L. Dubcova, P. Solin, J. Cerveny, P. Kus},
      title = {Space and Time Adaptive Two-Mesh hp-FEM for Transient Microwave Heating Problems},
      journal = {Electromagnetics},
      status = {accepted},
      year = {2009}
    }

::

    @Article{Hermes-fluid-mechanics,
      author = {P. Solin, J. Cerveny, L. Dubcova, I. Dolezel},
      title = {Multi-Mesh hp-FEM for Thermally Conductive Incompressible Flow},
      journal = {Proceedings of ECCOMAS Conference COUPLED PROBLEMS 2007 (M. Papadrakakis, E. Onate, 
                 B. Schrefler Eds.), CIMNE, Barcelona},
      year = {2007},
      pages = {677 - 680}
    }

Other papers that may be even closer to what you do can be found in the 
`publications section  <http://hpfem.org/publications/>`_ of the hp-FEM group home page.


