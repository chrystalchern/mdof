SHELL         = /bin/bash
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = docs
BUILDDIR      = _build

NOTEBOOKS = notebooks/00_Overview.ipynb \
            notebooks/01_SISO_Intro.ipynb \
            notebooks/02_SISO_Event.ipynb \
            notebooks/03_SISO_History.ipynb \
            notebooks/04_MIMO_Intro.ipynb \
	    	notebooks/06_MIMO_History.ipynb \


# Put this first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)


test:
	pytest --nbmake $(NOTEBOOKS)

publish: test
	mkdir -p site/
	cp -r _build/html/* site/
	touch site/.nojekyll
	git add site && git commit -m'cc - rebuild site'
	git subtree push --prefix site origin gh-pages
	make speak

speak:
	Say -v Whisper "dunzo, website is published"

html: Makefile
	mkdir -p ./docs/examples
	mkdir -p ./docs/examples/figures
	cp ./notebooks/0[0-9]_* ./docs/examples
	cp -r ./notebooks/figures/* ./docs/examples/figures
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
	sed -i .bak 's/\\bm{/\\boldsymbol{/g' $(BUILDDIR)/html/theory/srim.html
	

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

APIDOC = python3 tools/doc.py
APIDIR = docs/user/

conda:
	for i in 7 8 9 10; do conda mambabuild -c local -c conda-forge etc/conda --py 3.$i; done

pypa:
	sudo ./etc/pypa/docker-build
	

.PHONY: docs
