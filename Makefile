
all: real complex examples

real:
	@$(MAKE) -C src real
complex:
	@$(MAKE) -C src complex
release:
	@$(MAKE) -C src release
debug:
	@$(MAKE) -C src debug
profile:
	@$(MAKE) -C src profile
real-release:
	@$(MAKE) -C src real-release


examples:
	@$(MAKE) -C examples

clean:
	@$(MAKE) -C src clean
	rm -f lib/*
