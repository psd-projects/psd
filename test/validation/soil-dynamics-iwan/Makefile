SHELL := /usr/bin/zsh

preprocess:
	PSD_PreProcess -validation Iwan -activeplot

solve:
	PSD_Solve -np 1 PSD-validation-iwan.edp -v 0

run: preprocess solve

clean:
	@rm -rf *.edp *.data
