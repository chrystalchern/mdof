SHELL=/bin/bash

NOTEBOOKS = notebooks/01_SISO_Intro.ipynb \
            notebooks/02_SISO_Event.ipynb \
            notebooks/05_MIMO_Event.ipynb \
	    	notebooks/06_MIMO_History.ipynb \
			studies/PeakPicking.ipynb \
			studies/06_MIMO_History_All_CGS_Motions.ipynb

test:
	#pytest --nbmake notebooks/*.ipynb 
	pytest --nbmake $(NOTEBOOKS)




# notebook:
	# for notebook in ssid_demo.ipynb state_space_studies.ipynb impulse_response_studies.ipynb time_history_studies.ipynb transfer_function.ipynb; do \
# 	for notebook in ssid_demo.ipynb state_space_studies.ipynb; do \
# 	cp tests/$$notebook docs/examples; \
# 	cd docs/examples/ \
#         && cat <(printf "$$notebook\n==========================\n" ) <(pandoc $$notebook -t rst --title "SDOF"  --extract-media ./ --resource-path ../../tests/ ) > $${notebook/ipynb/rst} \
# 	&& cd - ; done
