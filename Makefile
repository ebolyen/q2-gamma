PLOT_SRC = q2_gamma/visualizers/plot/assets/src
PLOT_DST = q2_gamma/visualizers/plot/assets/dist

all: viz-plot

lint:
	q2lint
	flake8
	npm --prefix $(PLOT_SRC) run lint

test: all
	py.test

test-cov: all
	py.test --cov=q2_gamma


install: all
	python setup.py install

dev: all
	pip install -e .

clean: distclean
	rm -rf $(PLOT_SRC)/node_modules

distclean:
	rm -rf $(PLOT_DST)

.PHONY: viz-plot
viz-plot: $(PLOT_DST)

$(PLOT_SRC)/node_modules: $(PLOT_SRC)/package-lock.json
	npm --prefix $(PLOT_SRC) install

$(PLOT_DST): $(PLOT_SRC)/node_modules $(shell find $(PLOT_SRC)) | $(PLOT_DST)/licenses
	npm --prefix $(PLOT_SRC) run build

$(PLOT_DST)/licenses: $(shell find $(PLOT_SRC)/licenses)
	mkdir -p $(PLOT_DST)/licenses
	cp $(PLOT_SRC)/licenses/* $(PLOT_DST)/licenses/
