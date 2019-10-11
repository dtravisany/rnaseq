# rnaseq

Tutorial RNA-Seq

## Objetivos:

El siguiente tutorial tiene como objetivo introducir al alumno en el análisis bioinformático de 
experimentos de secuenciación provenientes de RNA-Seq.

Se realizarán las siguientes tareas:

- Control de Calidad de los reads.
- Alineamiento de los reads a un genoma de referencia.
- Conversión de los archivos a conteo por genes.
- Analisis de los conteos con DeSeq2.
- Anotación de nuestro experimento utilizando BioMart.
 
## Materiales:

Se utilizaran los datos de una Nota Técnica de BMC: [The bench scientist's guide to statistical analysisof RNA-Seq data](https://bmcresnotes.biomedcentral.com/track/pdf/10.1186/1756-0500-5-506).
 
Estos corresponden a transcriptoma de la soya creciendo en condiciones ambientales normales o de alto contenido de Ozono.

### Genoma de Referencia:

Como vimos en clase se necesitan dos datos:
El genoma de referencia y la anotación en un formato como el [GFF](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)
Se utilizará el genoma de referencia y anotación de la soya correspondiente a la versión [Wm82.a2.v1](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/515/GCF_000004515.5_Glycine_max_v2.1).

### Reads:

Los reads están alojados en las carpetas de su grupo en el servidor del curso.
De todas maneras, si los quiere descargar puede encontrar los links [acá](http://www.ncbi.nlm.nih.gov/sra/?term=SRP009826) y descargarlos utilizando la suite de NCBI-SRA [sra-tools](https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/).

### Software:

- Para el filtrado de los reads utilizaremos [TrimGalore](https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md).
- Para el alineamiento de los reads filtrados utilizaremos [STAR](https://github.com/alexdobin/STAR).
- Para generar el conteo a partir de los archivos BAM (que es la forma binaria de un [SAM](https://samtools.github.io/hts-specs/SAMv1.pdf)) 
utilizaremos [HTSeq](http://www-huber.embl.de/HTSeq/doc/overview.html).
- Para el analisis de expresión diferencial utilizaremos [DESeq2](http://bioconductor.org/packages/release/bioc/html/DESeq2.html).
- Para finalizar visualizaremos nuestros resultados en el [Integrative Genome Viewer](https://software.broadinstitute.org/software/igv/).
- Para finalizar visualizaremos nuestros resultados en [Artemis](https://software.broadinstitute.org/software/igv/).

# Inicio del Práctico

## Quality Check de los Reads

Debe realizar un quality Check de los reads utilizando el siguiente comando:

			trim-galore


Ese comando filtrara todos los reads que tengan un largo menor a 50bp y un [phred quality value](https://www.illumina.com/documents/products/technotes/technote_Q-Scores.pdf) mayor a 35.

Recordemos la tabla de phred Quality:

| Valor Phred| Probabilidad de base errónea | Precisión |  
| ----- | ---- |----|
| 10 | 1 en 10 | 90% |
| 20 | 1 en 100 | 99% |
| 30 | 1 en 1000 | 99.9% |
| 40 | 1 en 10000 | 99.99% |
| 50 | 1 en 100000 | 99.999% |





