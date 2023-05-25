SHELL=/bin/bash


notebook:
	# for notebook in ssid_demo.ipynb state_space_studies.ipynb impulse_response_studies.ipynb time_history_studies.ipynb transfer_function.ipynb; do \
	for notebook in ssid_demo.ipynb state_space_studies.ipynb; do \
	cp tests/$$notebook docs/examples; \
	cd docs/examples/ \
        && cat <(printf "$$notebook\n==========================\n" ) <(pandoc $$notebook -t rst --title "SDOF"  --extract-media ./ --resource-path ../../tests/ ) > $${notebook/ipynb/rst} \
	&& cd - ; done