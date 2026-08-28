SHELL := /usr/bin/zsh

preprocess-von-mises:
	./../../../src/psd-preprocessor/PSD_PreProcess -problem elasto_plastic -model von_mises -dimension 2 \
	-tractionconditions 1 -dirichletconditions 2 -useMfront && \
	find . -name ControlParameters.edp -exec sed -i -e "s#Dbc0On 2#Dbc0On 1#g" {} \; && \
	find . -name ControlParameters.edp -exec sed -i -e "s#macro  Dbc0Ux#//macro  Dbc0Ux 0.#g" {} \; && \
	find . -name ControlParameters.edp -exec sed -i -e "s#Dbc1On 2#Dbc1On 3#g" {} \; && \
	find . -name ControlParameters.edp -exec sed -i -e "s#macro  Dbc1Uy#//macro  Dbc1Uy#g" {} \; && \
	find . -name OtherParameters.edp -exec sed -i -e "s# real nRes, nRes0;# real nRes, nRes0;\n if(mpirank==0)\n  ofstream rr(\"result.data\");#g" {} \; && \
	find . -name LinearFormBuilderAndSolver.edp -exec sed -i -e "s#  //------------------Screen output ----------------------//#  //------------------Screen output ----------------------//\n  if(mpirank==0) rr << u(Ri,0) << \"\\\t\" << tl <<endl;#g" {} \; && \

preprocess-drucker-prager:
	./../../../src/psd-preprocessor/PSD_PreProcess -problem elasto_plastic -model drucker_prager -dimension 2 \
	-lagrange 2 -dirichletconditions 4  && \
	find . -name OtherParameters.edp -exec sed -i -e "s# real nRes, nRes0;# real nRes, nRes0;\n if(mpirank==0)\n  ofstream rr(\"result-drucker-prager.data\");#g" {} \; && \
	find . -name LinearFormBuilderAndSolver.edp -exec sed -i -e "s#  // --- screen output ---- //#  // --- screen output ---- //\n  if(mpirank==0) rr << i+1 << \"\\\t\" << normalizedPressure <<endl;#g" {} \; && \

preprocess-drucker-prager-mfront:
	./../../../src/psd-preprocessor/PSD_PreProcess -problem elasto_plastic -model drucker_prager -dimension 2 \
	-lagrange 2 -dirichletconditions 4 -useMfront && \
	find . -name OtherParameters.edp -exec sed -i -e "s# real nRes, nRes0;# real nRes, nRes0;\n if(mpirank==0)\n  ofstream rr(\"result-drucker-prager-mfront.data\");#g" {} \; && \
	find . -name LinearFormBuilderAndSolver.edp -exec sed -i -e "s#  // --- screen output ---- //#  // --- screen output ---- //\n  if(mpirank==0) rr << i+1 << \"\\\t\" << normalizedPressure <<endl;#g" {} \; && \

solve-von-mises:
	./../../../src/psd-solve/PSD_Solve -np 1 Main.edp -mesh ./../../../data/meshes/2D/Geo-Files/Gmsh/quater_cylinder.msh -v 0

solve-drucker-prager:
	./../../../src/psd-solve/PSD_Solve -np 1 Main.edp -mesh ./../../../data/meshes/2D/Geo-Files/Gmsh/strip_footing_geomechanics.msh -v 0 -ksp_rtol 1e-15

solve-drucker-prager-mfront:
	./../../../src/psd-solve/PSD_Solve -np 1 Main.edp -mesh ./../../../data/meshes/2D/Geo-Files/Gmsh/strip_footing_geomechanics.msh -v 0 -ksp_rtol 1e-15

gnuplot-von-mises:
	gnuplot gnuplot-compare.gnu -p	

gnuplot-drucker-prager:
	gnuplot gnuplot-compare-drucker-prager.gnu -p

gnuplot-drucker-prager-mfront:
	gnuplot gnuplot-compare-drucker-prager-mfront.gnu -p

run-von-mises: preprocess-von-mises solve-von-mises gnuplot-von-mises

run-drucker-prager: preprocess-drucker-prager solve-drucker-prager gnuplot-drucker-prager

run-drucker-prager-mfront: preprocess-drucker-prager-mfront solve-drucker-prager-mfront gnuplot-drucker-prager

clean:
	@rm -rf *.edp *.data *.csv *.eps
