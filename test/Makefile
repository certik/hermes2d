all:
	@$(MAKE) -C 01-test-vis
	@$(MAKE) -C 02-laplace
	@$(MAKE) -C 03-dirichlet
	@$(MAKE) -C 04-neumann
	@$(MAKE) -C 05-newton
#	@$(MAKE) -C 06-elasticity
	@$(MAKE) -C 07-hcurl-simple
	@$(MAKE) -C 08-hcurl-bessel
	@$(MAKE) -C 09-navier-stokes
	@$(MAKE) -C 10-heat-cond
	@$(MAKE) -C 11-multi-mesh-1
	@$(MAKE) -C 12-multi-mesh-2
	@$(MAKE) -C 13-h1-hang
#	@$(MAKE) -C 14-hcurl-hang
	@$(MAKE) -C 15-curved
	@$(MAKE) -C 17-adapt


clean:
	rm -f *.o *.bak core core.[0-9]*
	@$(MAKE) clean -C 01-test-vis
	@$(MAKE) clean -C 02-laplace
	@$(MAKE) clean -C 03-dirichlet
	@$(MAKE) clean -C 04-neumann
#	@$(MAKE) clean -C 06-elasticity
	@$(MAKE) clean -C 05-newton
	@$(MAKE) clean -C 07-hcurl-simple
	@$(MAKE) clean -C 08-hcurl-bessel
	@$(MAKE) clean -C 09-navier-stokes
	@$(MAKE) clean -C 10-heat-cond
	@$(MAKE) clean -C 11-multi-mesh-1
	@$(MAKE) clean -C 12-multi-mesh-2
	@$(MAKE) clean -C 13-h1-hang
#	@$(MAKE) clean -C 14-hcurl-hang
	@$(MAKE) clean -C 15-curved
	@$(MAKE) clean -C 17-adapt


