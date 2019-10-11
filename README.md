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
El genoma de referencia y la anotación en un formato como el [GFF](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md). 
En este caso, para la anotación utilizaremos el formato [Gene Transfer Format GTF](http://mblab.wustl.edu/GTF22.html) que es una derivación del formato 'GFF'. Al igual que el 'GFF' está tabulado, 
pero contiene algunas convenciones que son especificas del atributo gene del 'GFF'.

Se utilizará el genoma de referencia y anotación de la soya correspondiente a la versión [Wm82.a2.v1](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/515/GCF_000004515.5_Glycine_max_v2.1).

### Reads:

Los reads están alojados en las carpetas de su grupo en el servidor del curso.
De todas maneras, si los quiere descargar puede encontrar los links [acá](http://www.ncbi.nlm.nih.gov/sra/?term=SRP009826) y descargarlos utilizando la suite de NCBI-SRA [sra-tools](https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/).

### Software:

- Para el filtrado de los reads utilizaremos [TrimGalore](https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md).
- Para el alineamiento de los reads filtrados utilizaremos [STAR](https://github.com/alexdobin/STAR).
- Para generar el conteo a partir de los archivos 'BAM' (que es la forma binaria de un [SAM](https://samtools.github.io/hts-specs/SAMv1.pdf)) 
utilizaremos [HTSeq](http://www-huber.embl.de/HTSeq/doc/overview.html).
- Para el analisis de expresión diferencial utilizaremos [DESeq2](http://bioconductor.org/packages/release/bioc/html/DESeq2.html).
- Para finalizar visualizaremos nuestros resultados en el [Integrative Genome Viewer](https://software.broadinstitute.org/software/igv/).
- Para finalizar visualizaremos nuestros resultados en [Artemis](https://software.broadinstitute.org/software/igv/).

# Inicio del Práctico

## Quality Check de los Reads

Utilizaremos Trim Galore para hacer un quality check y filtrar los reads de las muestras:

			trim_galore --length 50  --quality 35 <READSFILE.fastq>


Ese comando filtrara todos los reads que tengan un largo menor a '50bp' y un [phred quality value](https://www.illumina.com/documents/products/technotes/technote_Q-Scores.pdf) mayor a '35'.

Recordemos la tabla de phred Quality:

| Valor Phred| Probabilidad de base errónea | Precisión |  
| ----- | ---- |----|
| 10 | 1 en 10 | 90% |
| 20 | 1 en 100 | 99% |
| 30 | 1 en 1000 | 99.9% |
| 40 | 1 en 10000 | 99.99% |
| 50 | 1 en 100000 | 99.999% |

El resultado del comando nos entregará los reads filtrados con extensión `trimmed.fq` y un log con extensión `fastq_trimming_report.txt`

## Alineamiento / Mapeo de los reads al genoma de referencia.

La carpeta `/home/<GRUPO>/RNA_SEQ/Genome/` (donde `<GRUPO>` es el identificador de su grupo) contiene los archivos del genoma de la soya, el archivo fasta de
los cromosomas y la anotación de los genes en formato 'GTF'. Estos fueron  descargados desde el [ftp de NCBI](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/515/GCF_000004515.5_Glycine_max_v2.1)

Para mapear los lreads al genoma de referencia los alineadores deben generar un índice del genoma.

	STAR --runMode genomeGenerate --genomeDir /home/<GRUPO>/RNA_SEQ/Genome/ --genomeFastaFiles /home/<GRUPO>/RNA_SEQ/Genome/Gmax_a2_v1.fna --sjdbGTFfile /common/RNASeq_Workshop/Soybean/gmax_genome/Gmax_275_Wm82.a2.v1.gene_exons --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 100 --runThreadN 8



