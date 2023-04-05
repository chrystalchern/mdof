SHELL=/bin/bash


notebook:
	# for notebook in system_id_summary.ipynb system_id_studies.ipynb transfer_function.ipynb; do \
	for notebook in system_id_summary.ipynb system_id_studies.ipynb; do \
	cp tests/$$notebook docs/examples; \
	cd docs/examples/ \
        && cat <(printf "$$notebook\n==========================\n" ) <(pandoc $$notebook -t rst --title "SDOF"  --extract-media ./ --resource-path ../../tests/ ) > $${notebook/ipynb/rst} \
	&& cd - ; done

