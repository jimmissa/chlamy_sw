The GO.zip archive contains three files: BP.gmt, CC.gmt, and MF.gmt. These are three GMT files enabling enrichment analysis on GO annotations on Chlamydomonas reinhardtii using the newly released v6.1 reference genome. I created this because current GO enrichment analysis software does not have any default recognition of the new gene annotations of the v6.1 reference genome, which is unsurprising because it only recently came out. For example, if you go to g:Profiler (https://biit.cs.ut.ee/gprofiler/), in the Options section you are able to select Chlamydomonas reinhardtii. Unfortunately, inputting a list of v6.1 gene names into the Query box will return no results, because g:Profiler currently recognizes v5.6 gene names.

One thing you might observe is that orthologous genes between v5.6 and v6.1 have the same name, except the v6.1 name has the string '_4532' appended to it. You might think, then, that you can simply use Bash or Python to remove all the '_4532' from the end of each gene name and just input that into g:Profiler or another software that recognizes v5.6 gene names. However, from my testing at least one third of gene names from the v6.1 gene names don't translate to v5.6 names in this manner, so this isn't the best way to do things. Instead, software like g:Profiler allows you to prepare your own custom GMT files for GO enrichment analysis, and this feature was introduced so that you could do enrichment analysis on non-model organisms. Chlamydomonas reinhardtii is a model organism, but obviously this custom file preparation is also very helpful while these software are catching up to the naming conventions of more up-to-date reference assemblies.

To do enrichment analysis on Chlamydomonas reinhardtii using these GMT files, you can:
1. Download the Zip archive
2. Go to g:Profiler (https://biit.cs.ut.ee/gprofiler/) — you can also use a different software but I'll focus on this one
3. On the right-half of the page, click on the black bar that says 'Bring your data (Custom GMT)'
4. Upload the GO.zip archive
5. Copy and paste your v6.1 gene names into the Query box on the left-half of the page, and finally hit 'Run query' to see the results




______________________________

The rest of this just explains how I prepared these GMT files. 

1) 

I first went to the Download options for the Chlamydomonas reinhardtii v6.1 genome on Phytozome (https://data.jgi.doe.gov/refine-download/phytozome?organism=CreinhardtiiCC-4532&expanded=707) and I downloaded the annotation info file ('CreinhardtiiCC_4532_707_v6.1.annotation_info.txt'). This file contains GO annotations of each gene in the 10th field. (Some genes have several GO annotations, some have none etc).

2)

I needed to also prepare a list of all possible GO annotations of each ontology type (BP = biological process, CC = cellular component, MF = molecular function). To do this, I used an R library called GO.db and I created four files: one file containing all GO annotations in the GO.db library and another three containing GO annotations only of each respective ontology. In the resulting files, the first field of each row was a GO annotation and the second tab-separated field was the corresponding vocabulary of that GO annotation. (This already brings us close to the GMT file format.) This is the script that will prepare you the same set of files:

suppressMessages(library(GO.db))

columns(GO.db)

go <- keys(GO.db, keytype="GOID")

df <- select(GO.db, columns=c("GOID","TERM"), keys=go, keytype="GOID")
df_bp <- select(GO.db, columns=c("GOID","TERM"), keys="BP", keytype="ONTOLOGY")
df_cc <- select(GO.db, columns=c("GOID","TERM"), keys="MF", keytype="ONTOLOGY")
df_mf <- select(GO.db, columns=c("GOID","TERM"), keys="CC", keytype="ONTOLOGY")

write.table(df, "df_goterms.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(df_bp[,2:3], "df_bp_goterms.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(df_cc[,2:3], "df_cc_goterms.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(df_mf[,2:3], "df_mf_goterms.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)

To give you an idea of how these files look like, a few lines of the df_bp_goterms.txt file look like this:

GO:0000379	tRNA-type intron splice site recognition and cleavage
GO:0000380	alternative mRNA splicing, via spliceosome
GO:0000381	regulation of alternative mRNA splicing, via spliceosome
GO:0000387	spliceosomal snRNP assembly

3) 

Before proceeding, I had to be sure that all the GO IDs in the v6.1 C. reinhardtii annotation info file were a subset of this supposedly complete set of all GO IDs as stored in the GO.db library. I extracted all GO IDs from the v6.1 annotation info file, and I also extracted the GO IDs of the df_goterms.txt file (which contains all GO annotations of all ontologies), made sure both were sorted, and then compared them using the comm command on Bash (with the -23 flags). As it happens, there are 50 or so GO annotations in the annotation info file that were not in the GO.db library (although this is out of like 50,000 or something in the GO.db library). I manually searched each one into geneontology.org and found that all of these annotations were either 1) obsolete or 2) had been replaced with a new GO term (which explains why they were not in the GO.db library). I made the decision to include these GO terms. If the GO annotation was obsolete, I manually added it to either df_bp_goterms.txt, df_cc_goterms.txt, or df_mf_goterms.txt depending on its ontology and I gave it its old vocabulary as specified on geneontology.org. The vocabulary includes the word 'obsolete' at the beginning so, if you're worried about this, you can just ignore any enriched vocabulary terms which begin with the word 'obsolete'. As for the GO IDs replaced with new GO IDs, I simply kept the old GO ID and added on the new vocabulary. Again, if you're worried about this, I reprint below the full list of replaced GO IDs so you can avoid them if you want to:

BP
GO:0006333	chromatin organization
GO:0006348	subtelomeric heterochromatin formation
GO:0006461	protein-containing complex assembly
GO:0006464	protein modification process
GO:0006827	iron ion transmembrane transport
GO:0007016	obsolete cytoskeletal anchoring at plasma membrane
GO:0007050	regulation of cell cycle
GO:0007067	mitotic cell cycle
GO:0009108	obsolete coenzyme biosynthetic process
GO:0009405	obsolete pathogenesis
GO:0015672	inorganic cation transmembrane transport
GO:0015696	ammonium transmembrane transport
GO:0015991	proton transmembrane transport
GO:0015992	proton transmembrane transport
GO:0016568	chromatin organization
GO:0035058	non-motile cilium assembly
GO:0042384	cilium assembly
GO:0043044	chromatin remodeling
GO:0051188	obsolete cofactor biosynthetic process
GO:0051297	centrosome cycle

MF
GO:0001104	transcription coregulator activity
GO:0003840	obsolete gamma-glutamyltransferase activity
GO:0004012	ATPase-coupled intramembrane lipid transporter activity
GO:0004584	obsolete dolichyl-phosphate-mannose-glycolipid alpha-mannosyltransferase activity
GO:0004871	obsolete signal transducer activity
GO:0004872	signaling receptor activity
GO:0005086	guanyl-nucleotide exchange factor activity
GO:0008144	obsolete drug binding
GO:0008536	small GTPase binding
GO:0008565	obsolete protein transporter activity
GO:0008601	protein phosphatase regulator activity
GO:0009540	zeaxanthin epoxidase activity
GO:0015002	obsolete heme-copper terminal oxidase activity
GO:0015238	xenobiotic transmembrane transporter activity
GO:0016820	ATPase-coupled transmembrane transporter activity
GO:0016876	aminoacyl-tRNA ligase activity
GO:0022891	transmembrane transporter activity
GO:0042623	ATP hydrolysis activity
GO:0043140	3'-5' DNA helicase activity
GO:0043141	5'-3' DNA helicase activity
GO:0048037	obsolete cofactor binding
GO:0050662	obsolete coenzyme binding
GO:0070011	peptidase activity
GO:0098519	obsolete nucleotide phosphatase activity, acting on free nucleotides

CC
GO:0030529	ribonucleoprotein complex
GO:0031513	non-motile cilium
GO:0043234	protein-containing complex


4) 

OK, so now I have 1) the annotation info file containing all the GO annotations of each v6.1 gene and 2) three files containing the GO IDs and vocabularies of each of the three ontologies. What I did next was write a Bash script that would take these files and spit out a GMT-formatted file, which effectively tells you on each row all the gene IDs which can be annotated by a specific given GO ID. I've uploaded this Bash script into the same folder as this README file, it's just called script.sh. To use the script, it has to be in the same directory which contains both the annotation file (which in my directory is named Cr.annotation_info.txt) and the df_bp/cc/mf_goterms.txt files. You can then use the script like this for each of the df_bp/cc/mf_goterms.txt files to create each of the three GMT files:

script.sh df_bp_goterms.txt > BP.gmt
script.sh df_cc_goterms.txt > CC.gmt
script.sh df_mf_goterms.txt > MF.gmt


I wrote this complete explanation based off of what I did from memory. If something doesn't make sense or doesn't work if you're trying to reproduce this yourself or change something, feel free to contact me.
