# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0002376","immune system process",0.5346209541723028,3.940548120738272,1,-0,"immune system process"),
c("GO:0002831","regulation of response to biotic stimulus",0.13539125829869797,6.603642527553255,0.7199218144362337,-0,"regulation of response to biotic stimulus"),
c("GO:0001952","regulation of cell-matrix adhesion",0.013778462205944099,3.561315969948284,0.7940447781136837,0.15173037,"regulation of response to biotic stimulus"),
c("GO:0002682","regulation of immune system process",0.2946264030637981,2.7230840431792416,0.8231267505663975,0.18910896,"regulation of response to biotic stimulus"),
c("GO:0007171","activation of transmembrane receptor protein tyrosine kinase activity",0.0005484791951702717,2.272969590019736,0.6758006962007334,0.46953855,"regulation of response to biotic stimulus"),
c("GO:0008285","negative regulation of cell population proliferation",0.09303204384988117,2.0751150310465656,0.73845444607297,0.6673656,"regulation of response to biotic stimulus"),
c("GO:0010511","regulation of phosphatidylinositol biosynthetic process",0.0009440490389597404,3.162175468643644,0.8442105479862082,0.12935154,"regulation of response to biotic stimulus"),
c("GO:0010543","regulation of platelet activation",0.005129111503925632,2.07696654849822,0.7736788983760255,0.53602191,"regulation of response to biotic stimulus"),
c("GO:0010628","positive regulation of gene expression",0.22275900791482073,2.035067769028446,0.7300974132785869,0.54404134,"regulation of response to biotic stimulus"),
c("GO:0021936","regulation of cerebellar granule cell precursor proliferation",0.0015789552588235097,2.3250810979952155,0.8493700674533035,0.13311865,"regulation of response to biotic stimulus"),
c("GO:0030155","regulation of cell adhesion",0.10961273612569519,2.6897234060933513,0.8275469242379686,0.17517661,"regulation of response to biotic stimulus"),
c("GO:0032101","regulation of response to external stimulus",0.3083018176633769,6.406430161402929,0.7058223978435804,0.68784959,"regulation of response to biotic stimulus"),
c("GO:0032103","positive regulation of response to external stimulus",0.11609143886137316,4.044745954885663,0.6415366313347473,0.63757584,"regulation of response to biotic stimulus"),
c("GO:0032648","regulation of interferon-beta production",0.007346297098947276,3.9034106337865286,0.7526601479941756,0.48582914,"regulation of response to biotic stimulus"),
c("GO:0032700","negative regulation of interleukin-17 production",0.002017738614959727,2.357374833613546,0.7164012136038872,0.6784416,"regulation of response to biotic stimulus"),
c("GO:0033120","positive regulation of RNA splicing",0.005780638305461228,2.33906933605301,0.7680895077951684,0.38657865,"regulation of response to biotic stimulus"),
c("GO:0034154","toll-like receptor 7 signaling pathway",0.0006149615218575774,2.1296353231063914,0.5941667117711902,0.48658216,"regulation of response to biotic stimulus"),
c("GO:0035790","platelet-derived growth factor receptor-alpha signaling pathway",0.0010603931106625253,2.012172978498394,0.7769886073067497,0.63430975,"regulation of response to biotic stimulus"),
c("GO:0038127","ERBB signaling pathway",0.013874861579640694,2.136069902428561,0.7470821157190035,0.61439365,"regulation of response to biotic stimulus"),
c("GO:0045071","negative regulation of viral genome replication",0.00546152313736216,5.032591971528533,0.7688439924117797,0.33691342,"regulation of response to biotic stimulus"),
c("GO:0045653","negative regulation of megakaryocyte differentiation",0.0006914161975479789,5.74806262927508,0.7012010097152271,0.12717242,"regulation of response to biotic stimulus"),
c("GO:0045932","negative regulation of muscle contraction",0.002360122597399351,2.5406678094463677,0.7307859335724403,0.51539674,"regulation of response to biotic stimulus"),
c("GO:0048008","platelet-derived growth factor receptor signaling pathway",0.005607784256074232,2.064751849985709,0.7584815432699447,0.69089045,"regulation of response to biotic stimulus"),
c("GO:0050792","regulation of viral process",0.020692624181423886,2.281475019285434,0.8487326604342231,0.15582003,"regulation of response to biotic stimulus"),
c("GO:0050818","regulation of coagulation",0.016680415765844992,2.7056339177927855,0.7648861878159862,0.57072124,"regulation of response to biotic stimulus"),
c("GO:0051209","release of sequestered calcium ion into cytosol",0.020646086552742773,2.962532401400879,0.73597955245238,0.62372036,"regulation of response to biotic stimulus"),
c("GO:0051282","regulation of sequestering of calcium ion",0.02990707466028445,2.7564252825411417,0.8261416210229949,0.159719,"regulation of response to biotic stimulus"),
c("GO:0051897","positive regulation of protein kinase B signaling",0.01250864976621656,2.605096352734407,0.6909741222389679,0.68426887,"regulation of response to biotic stimulus"),
c("GO:0060700","regulation of ribonuclease activity",0.003676472665808003,2.086180762791448,0.8374797441845514,0.17055739,"regulation of response to biotic stimulus"),
c("GO:0061000","negative regulation of dendritic spine development",0.0007213332445572665,2.512170631277455,0.7881604586384703,0.55274195,"regulation of response to biotic stimulus"),
c("GO:0061041","regulation of wound healing",0.020184034382266,5.18230937386111,0.7135882874322003,0.56376881,"regulation of response to biotic stimulus"),
c("GO:0150118","negative regulation of cell-substrate junction organization",0.0017351887265386777,2.3796636700406566,0.7710744176544369,0.35183395,"regulation of response to biotic stimulus"),
c("GO:1900246","positive regulation of RIG-I signaling pathway",0.001209978345708963,2.0134454455957767,0.6706353266429425,0.63758222,"regulation of response to biotic stimulus"),
c("GO:1901184","regulation of ERBB signaling pathway",0.014669325383553995,2.3407535945191826,0.745404934404013,0.55211003,"regulation of response to biotic stimulus"),
c("GO:1902093","positive regulation of flagellated sperm motility",0.0014725835361238204,2.515849518008878,0.7746967235774194,0.13259502,"regulation of response to biotic stimulus"),
c("GO:1903034","regulation of response to wounding",0.024332531567553873,4.710634332045057,0.7231549306831511,0.66062802,"regulation of response to biotic stimulus"),
c("GO:1904783","positive regulation of NMDA glutamate receptor activity",0.00019612286372755168,2.1715820067956995,0.7326489920524103,0.48708725,"regulation of response to biotic stimulus"),
c("GO:2000106","regulation of leukocyte apoptotic process",0.010018886631776963,2.8453959840100698,0.8441337965241322,0.14867302,"regulation of response to biotic stimulus"),
c("GO:2000107","negative regulation of leukocyte apoptotic process",0.00628590398828475,2.964766385771546,0.7769996704280198,0.60303671,"regulation of response to biotic stimulus"),
c("GO:2000342","negative regulation of chemokine (C-X-C motif) ligand 2 production",0.0003689769131145464,2.5740104350124193,0.7323130469788901,0.62093385,"regulation of response to biotic stimulus"),
c("GO:0006335","DNA replication-dependent chromatin assembly",0.008995058800792455,13.203689066069026,0.9577669450350526,0,"DNA replication-dependent chromatin assembly"),
c("GO:0001731","formation of translation preinitiation complex",0.010215009495504515,2.0943885317675472,0.9411729307931986,0.30167928,"DNA replication-dependent chromatin assembly"),
c("GO:0006325","chromatin organization",0.8645860379867383,2.072656076096586,0.957093019229528,0.41954381,"DNA replication-dependent chromatin assembly"),
c("GO:0006334","nucleosome assembly",0.08985751275056234,12.000799068681786,0.9357078200479761,0.57676998,"DNA replication-dependent chromatin assembly"),
c("GO:0032200","telomere organization",0.1347397314971624,2.744400801200614,0.9621373077261768,0.29927009,"DNA replication-dependent chromatin assembly"),
c("GO:0071824","protein-DNA complex subunit organization",0.20097939769205936,5.3755246292532455,0.9483702387333092,0.6035527,"DNA replication-dependent chromatin assembly"),
c("GO:0006396","RNA processing",4.052084515125943,6.970287705532774,0.9813981875499013,0.00645926,"RNA processing"),
c("GO:0000209","protein polyubiquitination",0.10754181164938563,2.4377255220202865,0.9888705592015952,0.11605825,"RNA processing"),
c("GO:0007548","sex differentiation",0.04039798581154129,2.1104684427242786,0.9545703508303669,0.47888703,"RNA processing"),
c("GO:0010467","gene expression",16.693722203531454,2.1047361491815377,0.9904954947692869,0.23325086,"RNA processing"),
c("GO:0018279","protein N-linked glycosylation via asparagine",0.012020004665064863,2.0836596739583095,0.9851926508594188,0.27717869,"RNA processing"),
c("GO:0032606","type I interferon production",9.97234900309585E-06,2.6563832359391175,0.9450643676178678,0.65126578,"RNA processing"),
c("GO:0032608","interferon-beta production",9.97234900309585E-06,3.9034106337865286,0.9450643676178678,0.15432422,"RNA processing"),
c("GO:0042713","sperm ejaculation",0.00018615051472445585,3.1110621682080826,0.9398237429473112,0.24880649,"RNA processing"),
c("GO:0042819","vitamin B6 biosynthetic process",0.12152304495172601,2.314466012377549,0.9878531566273702,0.19448401,"RNA processing"),
c("GO:0071608","macrophage inflammatory protein-1 alpha production",3.3241163343652828E-06,2.098917146010349,0.9469682645944189,0.62980563,"RNA processing"),
c("GO:0010761","fibroblast migration",0.0023435020157275242,2.220572671895782,0.9978305855666493,0.00404876,"fibroblast migration"),
c("GO:0019079","viral genome replication",0.025459407004903702,3.049306014373033,1,-0,"viral genome replication"),
c("GO:0030219","megakaryocyte differentiation",0.003324116334365283,5.0355934637599225,0.9647805982223916,0.00412086,"megakaryocyte differentiation"),
c("GO:0002305","CD8-positive, gamma-delta intraepithelial T cell differentiation",0.00035235633144272,2.488580030979042,0.8790476830597684,0.59236202,"megakaryocyte differentiation"),
c("GO:0046548","retinal rod cell development",0.0013529153480866703,2.2877559973054256,0.9343545958081513,0.48016263,"megakaryocyte differentiation"),
c("GO:0034109","homotypic cell-cell adhesion",0.005893658260829647,3.6324188335961303,0.9845604915535662,0.00424471,"homotypic cell-cell adhesion"),
c("GO:0031589","cell-substrate adhesion",0.04902074358288483,2.014449949948126,0.9846299083810062,0.56714532,"homotypic cell-cell adhesion"),
c("GO:0044419","biological process involved in interspecies interaction between organisms",1.030878281729696,3.965773578139414,1,-0,"biological process involved in interspecies interaction between organisms"),
c("GO:0051607","defense response to virus",0.15502681348579372,10.205394684570665,0.8116299889228237,-0,"defense response to virus"),
c("GO:0002227","innate immune response in mucosa",0.0019246633575974988,2.089978978173402,0.8312366159893443,0.60860947,"defense response to virus"),
c("GO:0006950","response to stress",4.85393118141787,2.7978759733306298,0.8570436752048456,0.47380831,"defense response to virus"),
c("GO:0006952","defense response",1.037047841646278,4.766127603367186,0.8412265592402113,0.52156191,"defense response to virus"),
c("GO:0007204","positive regulation of cytosolic calcium ion concentration",0.015726394377882154,4.008396441117039,0.81739548787806,0.53317028,"defense response to virus"),
c("GO:0007596","blood coagulation",0.0625133317840735,4.696624776171087,0.6403003312174979,0.41576626,"defense response to virus"),
c("GO:0009605","response to external stimulus",1.994426587106823,4.420599177132886,0.8674490587706031,0.38112452,"defense response to virus"),
c("GO:0009607","response to biotic stimulus",0.8976609955136728,5.054582966993217,0.8756706951518122,0.29406069,"defense response to virus"),
c("GO:0009611","response to wounding",0.11925599761168888,2.5767027646777976,0.8637143213259871,0.43610872,"defense response to virus"),
c("GO:0032486","Rap protein signal transduction",0.008842149449411653,2.281914605481856,0.7735913244295561,0.20808579,"defense response to virus"),
c("GO:0034097","response to cytokine",0.1986225992109944,4.8615445406174524,0.8365470801289745,0.30072489,"defense response to virus"),
c("GO:0034340","response to type I interferon",0.003696417363814195,5.357917893945226,0.7715899254559238,0.59774285,"defense response to virus"),
c("GO:0035455","response to interferon-alpha",0.002197240897015452,4.088778875025995,0.8507681762314137,0.52737122,"defense response to virus"),
c("GO:0050817","coagulation",0.06288895692985678,4.740139397244205,0.9419683158541938,0.49452235,"defense response to virus"),
c("GO:0050878","regulation of body fluid levels",0.08865750675385646,3.205132509376973,0.7993166770089344,0.59717874,"defense response to virus"),
c("GO:0061844","antimicrobial humoral immune response mediated by antimicrobial peptide",0.005355151414662471,2.7359445425568603,0.8258891357074413,0.61133323,"defense response to virus"),
c("GO:0070106","interleukin-27-mediated signaling pathway",0.0008875390612755306,2.738405743662763,0.7357243556751705,0.55543207,"defense response to virus"),
c("GO:0071357","cellular response to type I interferon",0.0025894866244705553,5.681563292937597,0.7692497028943367,0.61730497,"defense response to virus"),
c("GO:0071881","adenylate cyclase-inhibiting adrenergic receptor signaling pathway",0.00021939167806810867,2.450494411718233,0.8111759963720074,0.16863925,"defense response to virus"),
c("GO:0061644","protein localization to CENP-A containing chromatin",0.00011634407170278491,8.57309535754689,0.9680010862628416,0.0035197,"protein localization to CENP-A containing chromatin"),
c("GO:0051235","maintenance of location",0.07821978146394948,2.2354467044975688,0.9847188185031213,0.13071338,"protein localization to CENP-A containing chromatin"),
c("GO:0051651","maintenance of location in cell",0.04376531565825332,2.927032811824397,0.961055373617491,0.29558997,"protein localization to CENP-A containing chromatin"),
c("GO:0071887","leukocyte apoptotic process",0.0052487796919627825,2.5665029395498364,0.9908756733398647,0.00421905,"leukocyte apoptotic process"),
c("GO:0042771","intrinsic apoptotic signaling pathway in response to DNA damage by p53 class mediator",0.004417750608371461,2.3212703388412255,0.7541082738119219,0.599249,"leukocyte apoptotic process"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="C:/Yuchen/WCM/Courses/2023Spring/CMPB5004/project/angsd_project/checkpoint/final_0425/final_revigo_treemap.pdf", width=16, height=9 ) # width and height are in inches

# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = "Revigo TreeMap",
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  position.legend = "none"
)

dev.off()

