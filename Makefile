SHELL=bash
CONDA_ENV = amma

CONDA = $(shell which conda)
ifeq ($(CONDA),)
	CONDA=${HOME}/miniconda3/bin/conda
endif

default: help

clean: ## cleanup the project
	@rm -rf _site
	@rm -rf .sass-cache
	@rm -rf .bundle
	@rm -rf vendor
.PHONY: clean

create-env: ## create conda environment
	if ${CONDA} env list | grep '^${CONDA_ENV}'; then \
	    ${CONDA} env update -f environment.yml; \
	else \
	    ${CONDA} env create -f environment.yml; \
	fi
.PHONY: create-env

ACTIVATE_ENV = source $(shell dirname $(dir $(CONDA)))/bin/activate ${CONDA_ENV}
install: clean ## install dependencies for website
	$(ACTIVATE_ENV) && \
		cd docs && \
		gem install bundler && \
		bundle install
.PHONY: install

serve: ## run a local server
	$(ACTIVATE_ENV) && \
		cd docs && \
		bundle exec jekyll serve
.PHONY: serve

launch-jupyter: ## launch Jupyter
	$(ACTIVATE_ENV) && \
		jupyter notebook
.PHONY: launch-jupyter

generate-reports: ## generate HTML reports from the Jupyter Notebooks
	$(ACTIVATE_ENV) && \
		jupyter nbconvert --to=html src/*.ipynb --output-dir docs/ && \
		jupyter nbconvert --to=html src/full/*.ipynb --output-dir docs/full/ && \
		jupyter nbconvert --to=html src/sex-driven-aging/*.ipynb --output-dir docs/sex-driven-aging/ && \
		jupyter nbconvert --to=html src/microbiota-driven/*.ipynb --output-dir docs/microbiota-driven/
.PHONY: generate-reports


help:
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'
.PHONY: help
