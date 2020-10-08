# Receta para el análisis de rnaseq

Tutorial RNA-Seq y Expresión Diferencial

## Objetivos:

El siguiente tutorial tiene como objetivo introducir al alumno en el análisis bioinformático de 
experimentos de secuenciación provenientes de RNA-Seq.

Se realizarán las siguientes tareas:

- Control de Calidad de los reads.
- Alineamiento de los reads a un genoma de referencia.
- Conversión de los archivos a conteo por genes.
- Analisis de los conteos con DeSeq2.
- Enriquecimiento de Categorías GO.
 
## Materiales:

Una buena fuente de Datos de Expresion es el [Gene Expression Atlas](https://www.ebi.ac.uk/gxa/home) de [EBI](https://www.ebi.ac.uk)

En este tutorial utilizaremos datos descargados desde Genbank, encontrados en este Atlas.
Corresponde a datos del artículo [Reversal of Sepsis‐Like Features of Neutrophils by Interleukin‐1 Blockade in Patients With Systemic‐Onset Juvenile Idiopathic Arthritis](https://onlinelibrary.wiley.com/doi/full/10.1002/art.40442)
Puede descargar la tabla con la descripción de los datos [aquí](https://www.ebi.ac.uk/gxa/experiments-content/E-GEOD-103170/resources/ExperimentDesignFile.RnaSeq/experiment-design)

Corresponde a 12 RNA-Seq de 6 individuos humanos. 
Será su responsabilidad darle contexto a los Datos de la investigación.

Los datos ya se encuentran en el servidor y no es necesario descargarlos nuevamente, se ponen los link a disposición para que pueda identificar de donde se obtienen los datos en las bases de datos internacionales.


### Genoma de Referencia:

Como vimos en clase se necesitan dos datos:
El genoma de referencia y la anotación en un formato como el [GFF](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md). 
En este caso, para la anotación utilizaremos el formato [Gene Transfer Format GTF](http://mblab.wustl.edu/GTF22.html) que es una derivación del formato 'GFF'. Al igual que el 'GFF' está tabulado, 
pero contiene algunas convenciones que son especificas del atributo gene del 'GFF'.

En el artículo, se utilizó la versión [GRCh37] del Genoma Humano, en nuestro caso utilizaremos la última versión del Genoma correspondiente a la [GRCh38.p13](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.39) liberado en febrero del 2019. El ftp con la data la puede encontrar [acá](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13)

### Reads:

Los reads están alojados en las carpetas de su grupo en el servidor del curso.
De todas maneras, si los quiere descargar puede encontrar los id en la [tabla](https://www.ebi.ac.uk/gxa/experiments-content/E-GEOD-103170/resources/ExperimentDesignFile.RnaSeq/experiment-design) \(SRR5984243-SRR5984254\) y descargarlos utilizando la suite de NCBI-SRA [sra-tools](https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/), especificamente la herramienta [fasterq-dump](https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump).

### Software:

- Para el filtrado de los reads utilizaremos [TrimGalore](https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md).
- Para el alineamiento de los reads filtrados utilizaremos [STAR](https://github.com/alexdobin/STAR).
- Para generar el conteo a partir de los archivos 'BAM' (que es la forma binaria de un [SAM](https://samtools.github.io/hts-specs/SAMv1.pdf)) 
utilizaremos [HTSeq](https://htseq.readthedocs.io/en/release_0.11.1/).
- Para el analisis de expresión diferencial utilizaremos [DESeq2](http://bioconductor.org/packages/release/bioc/html/DESeq2.html).
- Visualizaremos nuestros resultados en el [Integrative Genome Viewer](https://software.broadinstitute.org/software/igv/).
- Exploraremos un análisis de Enriquecimiento de las categorias [GO](http://geneontology.org/docs/ontology-documentation/)

# Inicio del Práctico

## Quality Check de los Reads

Utilizaremos Trim Galore para hacer un quality check y filtrar los reads de las muestras, entonces para cada grupo de reads ejecutaremos
la instrucción:

			trim_galore --length 50  --quality 35 <READSFILE.fastq>


Ese comando filtrara todos los reads que tengan un largo menor a `50bp` y un [phred quality value](https://www.illumina.com/documents/products/technotes/technote_Q-Scores.pdf) mayor a `35`.

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

### Generación de un índice del Genoma

Para mapear los reads al genoma de referencia los alineadores deben generar un índice del genoma. Esto permite que el proceso de mapear millones de reads se haga de forma eficiente.

La carpeta `/home/<GRUPO>/RNA_SEQ/Genome/` (donde `<GRUPO>` es el identificador de su grupo) contiene los archivos del genoma de la soya, el archivo fasta de
los cromosomas y la anotación de los genes en formato `GTF`. Estos fueron  descargados desde el [ftp de NCBI](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/515/GCF_000004515.5_Glycine_max_v2.1)


 `STAR` genera el índice  de la siguiente manera:

	STAR --runMode genomeGenerate --genomeDir /home/<GRUPO>/RNA_SEQ/Genome/ --genomeFastaFiles /home/<GRUPO>/RNA_SEQ/Genome/Gmax_a2_v1.fna --sjdbGTFfile /home/<GRUPO>/RNA_SEQ/Genome/Gmax_a2_v1.gtf --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 99 --runThreadN 8

Todos los resultados quedarán guardados en la dirección dada al parámetro `--genomeDir`

Con respecto a las opciones/parámetros utilizadas:

| Parámetro | Descripción |
| ---- | ---- |
| `--runMode` | indica el tipo de opcion que utilizará STAR, en este caso queremos generar un índice del genoma por lo que utilizamos la flag `genomeGenerate`|
| `--genomeDir` | indica donde se guardaran los resultados del indice y la ubicación de los archivos del genoma |
| `--genomeFastaFiles` | indica donde estan almacenadas las secuencias del genoma en formato `FASTA` |
| `--sjdbGTFfile` | sj: splice junction db: database GTFfile: archivo GTF indica la ubicación del archivo GTF para mejorar e improvisar el mapeo dado el modelo de los genes |
| `--sjdbOverhang` | Especifica el largo a considerar de la secuencia genómica alrededor del splice junction, este valor esta ligado al largo de los reads y deberia ser `max(ReadLength) - 1`  |
| `--runThreadN` | total de hebras que se ejecutaran en paralelo, este número no debe sobrepasar la cantidad de cores que tiene un computador y pruebas de escalamiento deberian ser ejecutadas para calcular el óptimo |


### Mapping

Dado que se genera el índice de los cromosomas y el genoma, ahora debemos proceder a mapear los reads que pasaron los filtros de calidad.

Para eso `STAR` utiliza:

	STAR --genomeDir /home/<GRUPO>/RNA_SEQ/Genome/ --readFilesIn READSFILE_sra.trimmed.fq --runThreadN 8 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix READSFILE_sra.trimmed

| Parámetro | Descripción |
| ---- | ---- |
| `--readFilesIn` | archivo de reads a mapear |
| `--genomeDir` | indica donde esta alojado el genoma |
| `--runThreadN`| cantidad de hebras |
| `--outSAMType`| tipo de archivo `SAM` o `BAM` en este caso indicamos `BAM`|
| `SortedByCoordinate`| Ordenar el archivo BAM por las coordenadas del genoma |
| `--outFileNamePrefix`| Prefijo para los archivos de salida |

Para mayor información sobre `STAR` puede visitar el [manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf) 


Este proceso generará un archivo `BAM` que tiene los reads mapeados al genoma.

##Generar los conteos

El `BAM` recien generado tiene las posiciones y el detalle de cada read mapeado a nuestro genoma de referencia,
pero para calcular las etadísticas del mapeo necesitamos un programa que nos permita transformar esta información detallada
en un simple archivo de abundancia de reads por gen. Para esto utilizaremos [HTSeq](https://htseq.readthedocs.io/en/release_0.11.1/).


		 htseq-count -s no -r pos -t exon -f bam  <BAMFILE> <GTF> > <SALIDA>

Ok, acabamos de terminar de trabajar con los grandes volumenes de datos de secuenciación y hemos terminado con un archivo 
pequeño que contiene los reads mapeados a cada genoma. Esto puede ser trabajado en su computador personal.

## Instalación de R en su computador

#### Nota: Si ya tiene instalado R, en su computador puede omitir este paso.

R es un lenguaje y ambiente para realizar estadística computacional. Puede obtener más información desde la [página del Proyecto](https://www.r-project.org/)



Descargar [R-base](https://www.r-project.org/)

Descargar [RStudio](https://www.rstudio.com/products/rstudio/#Desktop)

R tiene diversos webservers y servidores que guardan información y versiones de código que ya esta programado, este repositorio se llama [CRAN - Comprehensive R Archive Network](https://cran.r-project.org).

Para bioinformática existe un repositorio particular de código de R que se llama [bioconductor](https://www.bioconductor.org/)

## Cálculo de expresión diferencial

Para calcular la expresión diferencial necesitaremos la librería de R [DESeq2](https://www.bioconductor.org/packages/release/bioc/html/DESeq2.html). La publicación que detalla la metodología implementada en DESeq2 la puede encontrar [acá](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4302049/pdf/13059_2014_Article_550.pdf).

### Jupyter-Notebooks

He diseñado [Jupyter Notebooks](https://jupyter.org/) para hacer la clase más dinámica, generalmente los jupyter-notebook se cargan bien en Github, pero si tuviesen algun problema pueden copiar el link del Jupyter-Notebook en:

[https://nbviewer.jupyter.org/](https://nbviewer.jupyter.org/)













