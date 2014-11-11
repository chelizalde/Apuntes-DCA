TARGET=Matematicas
.PHONY: all book clean

all: book clean

book:
	@${MAKE} clean
	@echo "==============================================="
	@echo "Compilando con pdflatex, el archivo:" $(TARGET)
	@pdflatex -interaction=batchmode $(TARGET).tex
	@pdflatex -interaction=batchmode $(TARGET).tex
	@pdflatex -interaction=batchmode $(TARGET).tex

clean:
	@echo "==============================================="
	@echo "Quitando archivos irrelevantes."
	@rm -rf *.aux *.log *.fdb_latexmk *.synctex.gz *.out *.tdo *toc
