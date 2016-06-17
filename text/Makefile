word: ms.docx

ms.docx: ms.md
	pandoc --standalone --reference-docx=style.docx \
		--filter pandoc-citeproc --bibliography=bib/bib.bibtex \
		--csl=bib/g3.csl --from=markdown --to=docx --output=ms.docx ms.md