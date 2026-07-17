# MethylHome

## Table des matières

- [Description](#description)
- [Informations importantes](#informations-importantes)
- [Vue d'ensemble du pipeline](#vue-densemble-du-pipeline)
- [Paramètres d'entrée](#paramètres-dentrée)
  - [Paramètres obligatoires](#paramètres-obligatoires)
    - [`sample_sheet`](#sample_sheet)
  - [Paramètres optionnels](#paramètres-optionnels)
    - [`qc_ref_set`](#qc_ref_set)
    - [`CNV_focal`](#cnv_focal)
    - [`ref_m`](#ref_m)
    - [`ref_f`](#ref_f)
    - [`anno`](#anno)
    - [`file_data`](#file_data)
    - [`publish`](#publish)
    - [`output`](#output)
    - [`all_qc_metrics_file`](#all_qc_metrics_file)
    - [`all_qc_metrics_gs_file`](#all_qc_metrics_gs_file)
    - [`all_CNV_detail_file`](#all_cnv_detail_file)
    - [`all_CNV_metrics_file`](#all_cnv_metrics_file)
    - [`all_CNV_segment_file`](#all_cnv_segment_file)
    - [`all_tumor_purity_file`](#all_tumor_purity_file)
    - [`all_mgmt_file`](#all_mgmt_file)
- [Chargement des données](#chargement-des-données)
  - [`load_idats`](#load_idats)
    - [Entrée](#entrée)
    - [Sortie](#sortie)
  - [`load_idats_minfi`](#load_idats_minfi)
      - [Entrée](#entrée-1)
      - [Sortie](#sortie-1)
- [Module QC](#module-qc)
  - [`plot_qc`](#plot_qc)
      - [Entrée](#entrée-2)
      - [Sortie](#sortie-2)
  - [`merge_qc_metrics`](#merge_qc_metrics_gs)
      - [Entrée](#entrée-3)
      - [Sortie](#sortie-3)
  - [`compute_qc_gs`](#compute_qc_gs)
      - [Entrée](#entrée-4)
      - [Sortie](#sortie-4)
  - [`merge_qc_metrics_gs`](#merge_qc_metrics_gs)
      - [Entrée](#entrée-5)
      - [Sortie](#sortie-5)
  - [`plot_qc_gs`](#plot_qc_gs)
      - [Entrée](#entrée-6)
      - [Sortie](#sortie-6)
  - [`extract_xy_intensities`](#extract_qc_metrics)
      - [Entrée](#entrée-7)
      - [Sortie](#sortie-7)
  - [`merge_xy_intensities`](#merge_xy_intensities)
      - [Entrée](#entrée-8)
      - [Sortie intermédiaire](#sortie-intermédiaire)
  - [`control_sex`](#control_sex)
      - [Entrée](#entrée-9)
      - [Sortie](#sortie-8)
        - [`control_sex_report`](#control_sex_report)
        - [`all_predicted_sex.tsv`](#all_predicted_sextsv)
- [Module CNV](#module-cnv)
    - [Gestion des chromosomes sexuels](#gestion-des-chromosomes-sexuels)
    - [Particularité des puces EPICv2](#particularité-des-puces-epicv2)
    - [Jeu de données de référence](#jeu-de-données-de-référence)
    - [Recommandations](#recommandations)
      - [Construction d’un jeu de référence](#construction-dun-jeu-de-référence)
      - [Qualité des données](#qualité-des-données)
      - [Paramètres de binning](#paramètres-de-binning)
  - [`compute_cnv`](#compute_cnv)
    - [Entrée](#entrée-10)
    - [Sortie](#sortie-9)
      - [Visualisations](#visualisations)
      - [Fichiers de données](#fichiers-de-données)
  - [`merge_cnv_detail`](#merge_cnv_detail)
    - [Entrée](#entrée-11)
    - [Sortie](#sortie-10)
  - [`merge_cnv_segment`](#merge_cnv_segment)
    - [Entrée](#entrée-12)
    - [Sortie](#sortie-11)
  - [`merge_cnv_metrics`](#merge_cnv_metrics)
    - [Entrée](#entrée-13)
    - [Sortie](#sortie-12)
- [Module MGMT](#module-mgmt)
  - [Description](#description-1)
  - [`mgmt`](#mgmt)
    - [Entrée](#entrée-14)
    - [Sortie](#sortie-13)
  - [`merge_mgmt`](#merge_mgmt)
      - [Entrée](#entrée-15)
      - [Sortie](#sortie-14)
- [Module Tumor Purity](#module-tumor-purity)
  - [Description](#description-2)
  - [`compute_tumor_purity`](#compute_tumor_purity)
    - [Entrée](#entrée-16)
    - [Sortie](#sortie-15)
  - [`merge_tumor_purity`](#merge_tumor_purity)
      - [Entrée](#entrée-17)
      - [Sortie](#sortie-16)
- [Ressources supplémentaires](#ressources-supplémentaires)
  - [Où trouver les informations d'annotation dans les différents objets :](#où-trouver-les-informations-dannotation-dans-les-différents-objets-)
- [Preprocessing](#preprocessing)
  - [Normalisation avec `minfi`](#normalisation-avec-minfi)
  - [Normalisation avec `ewastools`](#normalisation-avec-ewastools)
    - [**Objets bruts** (*raw objects*)](#objets-bruts-raw-objects)
    - [**Beta-values**](#beta-values)


## Description

MethylHome est un pipeline Nextflow dédié à l’analyse de fichiers IDAT issus de puces de méthylation. Il permet d’extraire les métriques de contrôle qualité, de générer des rapports QC individuels, de réaliser un contrôle d’identitovigilance basé sur le sexe et de calculer les variations du nombre de copies (CNV).

## Informations importantes

L’inférence des CNV focaux nécessite un accès à Internet. Afin de permettre une exécution hors ligne, une option désactivant le calcul des CNV focaux ainsi que la génération des graphiques associés aux régions d’intérêt est disponible.

## Vue d’ensemble du pipeline

![workflow du pipeline MethylHome](MethylHome_workflow.png)
*Schéma global du workflow.*

## Paramètres d’entrée

### Paramètres obligatoires

Le pipeline nécessite deux paramètres d’entrée obligatoires.

#### `sample_sheet`

Fichier au format CSV ou TSV contenant au minimum les colonnes suivantes :

* **Sample_Name** : nom de l’échantillon. Cette valeur doit être unique.
* **Sample_IDAT** : combinaison du `Sentrix_ID` et du `Sentrix_Position` de l’échantillon (exemple : `207716530108_R04C01`).
* **Gender** : sexe de l’échantillon. Les valeurs autorisées sont `M`, `F` ou `U` lorsque l’information n’est pas disponible.
* **filepath** : chemin relatif ou absolu vers les fichiers IDAT, sans le suffixe `_Grn.idat` ou `_Red.idat`.

Exemple :

| Sample_Name | Sample_IDAT                    | Gender | file_path                              |
| ----------- | ------------------------------ | ------ | -------------------------------------- |
| GSM8461134  | GSM8461134_207716530108_R01C01 | M      | path/to/GSM8461134_207716530108_R01C01 |

Si la `sample_sheet` contient des lignes vides avant **Sample_Name**, le piepline ne fonctionnera pas proprement. Il est donc important d'insérer au minimum un caractère par ligne.


### Paramètres optionnels

Les paramètres suivants sont facultatifs ou disposent d’une valeur par défaut modifiable par l’utilisateur.

#### `qc_ref_set`

Fichier TSV contenant les valeurs de métriques QC d’un jeu de données de référence. Par défaut, il s’agit du jeu de données GEO [GSE274910](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274910). Cette référence est utilisée pour positionner les échantillons analysés au sein des distributions observées dans la cohorte de référence.

#### `CNV_focal`

*TRUE* ou *FALSE*, détermine si le calcul des CNV focaux et les représentations associées sont effectuées par le pipeline. Ces opérations nécessitant une connexion internet, ce paramètre doit être mis à la valeur FALSE pour fonctionner offline. Par défaut, les CNV focaux sont calculés : *TRUE*.  

#### `ref_m`

Objet `.Rdata` utilisé pour le calcul des CNV des échantillons masculins. Il contient la référence utilisée pour la normalisation et la réduction du bruit de fond des intensités des sondes. Par défaut, les échantillons non tumoraux du dépôt GEO [GSE306226](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE306226) sont utilisés.

#### `ref_f`

Objet `.Rdata` utilisé pour le calcul des CNV des échantillons féminins. Il contient la référence utilisée pour la normalisation et la réduction du bruit de fond des intensités des sondes. Par défaut, les échantillons non tumoraux du dépôt GEO [GSE306226](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE306226) sont utilisés.

#### `ref_mf`

Objet `.Rdata` utilisé pour le calcul des CNV des échantillons dont le sexe n'est pas renseigné. Il contient la référence utilisée pour la normalisation et la réduction du bruit de fond des intensités des sondes. Par défaut, les échantillons non tumoraux du dépôt GEO [GSE306226](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE306226) sont utilisés.

#### `anno`

Objet `.Rdata` contenant l’annotation (genomic binning) utilisée pour le calcul des CNV. Cette annotation permet la prise en charge des puces :

* Infinium Human Methylation 450K BeadChip ;
* Infinium MethylationEPIC v1.0 BeadChip ;
* Infinium MethylationEPIC v2.0 BeadChip.

Par défaut, l’annotation utilisée est de type EPIC/EPICv2.

#### `publish`

Mode de publication des résultats. Les valeurs possibles sont :

* `copy`
* `copyNoFollow`
* `link`
* `move`
* `rellink`
* `symlink`

La valeur par défaut est `copy`.

#### `output`

Répertoire de sortie dans lequel les résultats seront publiés. Par défaut ce dossier s'appelle : *output*.

#### `all_qc_metrics_file` 

Nom du fichier sans extension dans lequel vont être stockées les valeurs des métriques de contrôle qualité normalisées de tous les échantillons présents dans la `sample_sheet`. Par défaut : *all_qc_metrics*.


#### `all_qc_metrics_gs_file` 

Nom du fichier sans extension dans lequel vont être stockées les valeurs des métriques de contrôle qualité brutes de tous les échantillons présents dans la `sample_sheet`. Par défaut : *all_qc_metrics_gs*.

#### `all_CNV_detail_file`

Nom du fichier sans extension dans lequel vont être stockées les valeurs des régions d'intérêts des CNV de tous les échantillons présents dans la `sample_sheet`. Par défaut : *all_CNV_detail*.


#### `all_CNV_metrics_file`

Nom du fichier sans extension dans lequel vont être stockées les valeurs des métriques du calcul des CNV de tous les échantillons présents dans la `sample_sheet`. Par défaut : *all_CNV_metrics*.

#### `all_CNV_segment_file`

Nom du fichier sans extension dans lequel vont être stockées les valeurs des segments des CNV de tous les échantillons présents dans la `sample_sheet`. Par défaut : *all_CNV_segment*.

#### `all_tumor_purity_file`

Nom du fichier sans extension dans lequel vont être stockées les valeurs estimées de la pureté tumorale de tous les échantillons présents dans la `sample_sheet`. Par défaut : *all_tumor_purity*.

#### `all_mgmt_file`

Nom du fichier sans extension dans lequel vont être stockées les valeurs et le statut de méthylation du promoteur du gène *MGMT* de tous les échantillons présents dans la `sample_sheet`. Par défaut : *all_mgmt*.

## Chargement des données

### `load_idats`

Ce process lit les fichiers IDAT et génère un fichier `.rds` à partir de l’objet créé par la fonction `read_idats` du package `ewastools`.

#### Entrée

Channel Nextflow contenant, pour chaque échantillon, le tuple suivant :

* **sample_id** : identifiant de l’échantillon ;
* **sex** : sexe de l'échantillon ; 
* **file_green** : fichier IDAT vert ;
* **file_red** : fichier IDAT rouge.

#### Sortie

Tuple contenant :

* **sample_id** ;
* **sex** ;
* fichier `.rds` généré nommé `sample_id.rds`.

### `load_idats_minfi`

Ce process lit les fichiers IDAT et génère un fichier `.rds` à partir de l’objet créé par la fonction `read_metharray` du package `minfi`.

#### Entrée

Channel Nextflow contenant, pour chaque échantillon, le tuple suivant :

* **sample_id** : identifiant de l’échantillon ;
* **sex** : sexe de l'échantillon ; 
* **file_green** : fichier IDAT vert ;
* **file_red** : fichier IDAT rouge.

#### Sortie

Tuple contenant :

* **sample_id** ;
* **sex** ;
* fichier `.rds` généré nommé `sample_id_minfi.rds`.


## Module QC

![workflow du module QC du pipeline MethylHome](QC_workflow.png)
*Schéma du workflow du module QC.*

Le contrôle qualité des échantillons repose sur l’analyse combinée de plusieurs indicateurs techniques calculés à partir des fichiers bruts générés par les scanners Illumina iScan ou Illumina NextSeq 550.

Les métriques sont calculées à l’aide du package R `ewastools`, directement à partir des données natives produites par les puces Illumina.

Le module comprend :

* l’extraction et la représentation des valeurs des sondes de contrôle qualité telles qu’elles sont utilisées par GenomeStudio ;
* le calcul des métriques définies par le logiciel BeadArray Controls Reporter ;
* l’évaluation des intensités logarithmiques des signaux méthylés et non méthylés ;
* l’analyse de la fonction de distribution cumulative empirique (ECDF) des p-values de détection ;
* l’étude de la distribution des beta-values selon le type de sonde ;
* la comparaison des résultats à une cohorte de référence ;
* la génération d’un tableau récapitulatif des métriques, des seuils recommandés et de leur statut.

À l’exception de la représentation des beta-values, aucune normalisation n’est appliquée avant l’évaluation des métriques de contrôle qualité. Ce choix est motivé par les éléments suivants :

1. Les métriques de contrôle qualité définies par Illumina sont calculées sur les données brutes.
2. L’absence de normalisation évite qu’un échantillon borderline soit artificiellement amélioré et considéré comme conforme.

Les seuils d’interprétation utilisés correspondent aux recommandations officielles d’Illumina.

L’ensemble des valeurs calculées est exporté dans un fichier TSV afin de permettre leur réutilisation en dehors du module QC.


### `plot_qc`

Ce process calcule les métriques de contrôle qualité de chaque échantillon.

Deux catégories de métriques sont produites :

* **QC globaux** :

  * intensités log2 méthylées ;
  * intensités log2 non méthylées ;
  * taux de détection.

* **QC wet lab** :

  * 21 métriques de contrôle qualité décrites dans le guide *BeadArray Controls Reporter Software Guide* d’Illumina.

Il génère également un rapport PDF de contrôle qualité pour chaque échantillon.

#### Entrée

Tuple issu du process `load_idats` :

* **sample_id** ;
* **sex** ;
* **meth_rds**.

Le process utilise également :

* **qc_ref_set** : base de référence contenant les distributions des métriques QC fournie en entrée du pipeline.

#### Sortie

Par échantillon : 
* fichier TSV `sample_id_qc_metrics_output.tsv` contenant l’ensemble des métriques calculées.

* Rapport PDF comprenant :

  1. **Métriques QC BeadArray Illumina**

    Représentation des 21 métriques sous forme de distributions de densité calculées à partir du jeu de référence. La valeur de l’échantillon est indiquée par un point rouge. Chaque métrique est normalisée par sa propre variance.

  2. **Log-intensités méthylées et non méthylées**

    Visualisation des intensités logarithmiques des signaux méthylés (M) et non méthylés (U), permettant d’identifier une baisse globale d’intensité ou un déséquilibre entre les canaux.

  3. **ECDF des p-values de détection**

    Représentation de la fonction de distribution cumulative empirique des p-values de détection. Une proportion élevée de sondes présentant des p-values importantes peut traduire une qualité insuffisante de l’échantillon.

  4. **Distribution des beta-values**

    Histogrammes des beta-values avec densités séparées pour les sondes de type I et de type II. Cette représentation permet de vérifier la bimodalité attendue et de détecter d’éventuelles anomalies globales du profil de méthylation.

Le rapport contient également un tableau récapitulatif présentant les valeurs calculées ainsi que leur statut de validation.

Lorsque la qualité d’un échantillon est insuffisante (`colSums(is.na(meth_QC[["M"]] + meth_QC[["U"]])) > 1000`), certaines figures ne sont pas générées, notamment la distribution des beta-values et le taux de détection. Un avertissement est alors affiché au début du rapport PDF.

### `merge_qc_metrics`

Ce process concatène les fichiers produits par `plot_qc` en un unique fichier TSV.

Les lignes du fichier final sont triées selon le nom des échantillons.

#### Entrée

* all_qc_metrics ;
* Ensemble des fichiers TSV produits par `plot_qc`.

#### Sortie

Fichier :

* `all_qc_metrics.tsv`

contenant les métriques de contrôle qualité de tous les échantillons présents dans la `sample_sheet`.

### `compute_qc_gs`

Ce process extrait les valeurs brutes des sondes de contrôle qualité selon la même approche que celle utilisée par le logiciel GenomeStudio d’Illumina.

#### Entrée

Tuple issu du process `load_idats` :

* **sample_id** ;
* **sex** ;
* **meth_rds**.

#### Sortie

Pour chaque échantillon :

* fichier TSV `sample_id_qc_metrics_gs.tsv` contenant l’ensemble des valeurs brutes extraites.

### `merge_qc_metrics_gs`

Ce process concatène les fichiers générés par `compute_qc_gs` en un unique fichier TSV.

Les lignes du fichier final sont triées selon le nom des échantillons.

#### Entrée

* all_qc_metrics_gs ;
* Ensemble des fichiers TSV produits par `compute_qc_gs`.

#### Sortie

Fichier :

* `all_qc_metrics_gs.tsv`

contenant les valeurs brutes des sondes de contrôle qualité pour l’ensemble des échantillons présents dans la `sample_sheet`.

### `plot_qc_gs`

Ce process génère un rapport PDF récapitulant les valeurs brutes des sondes de contrôle qualité pour tous les échantillons.

#### Entrée

* **merged_qc_gs** : fichier TSV issu du process `merge_qc_metrics_gs`, contenant les valeurs brutes de contrôle qualité de l’ensemble des échantillons de la `sample_sheet`.

#### Sortie

Fichier PDF :

* `genome_studio_like_plot.pdf`

Ce rapport contient les représentations des métriques suivantes :

* Restoration ;
* Bisulfite Conversion I Green ;
* Bisulfite Conversion I Red ;
* Bisulfite Conversion II Green ;
* Bisulfite Conversion II Red ;
* Extension Green ;
* Extension Red ;
* Hybridization Green ;
* Hybridization Red ;
* Non-polymorphic Green ;
* Non-polymorphic Red ;
* Specificity I Green ;
* Specificity I Red ;
* Specificity II Green ;
* Specificity II Red ;
* Staining Green ;
* Staining Red ;
* Target Removal Green ;
* Target Removal Red.

À l’exception de la métrique *Restoration*, une figure distincte est générée pour chacun des canaux vert (*Green*) et rouge (*Red*).

### `extract_xy_intensities`

Ce process extrait les intensités des chromosomes X et Y pour chaque échantillon.

#### Entrée

Tuple issu du process `load_idats` :

* **meth_rds** ;
* **sex** ;
* **sample_id**.

#### Sortie

Pour chaque échantillon :

* fichier TSV `sample_id_xy_intensities.tsv` contenant les intensités des chromosomes X et Y ainsi que le sexe renseigné dans la `sample_sheet`.

Lorsque la qualité d’un échantillon est insuffisante (`colSums(is.na(meth_QC[["M"]] + meth_QC[["U"]])) > 1000`), les intensités X et Y ne sont pas prédites et le sexe est fixé à `"U"`.

### `merge_xy_intensities`

Ce process concatène les fichiers générés par `extract_xy_intensities` en un unique fichier TSV.

Les lignes du fichier final sont triées selon le nom des échantillons.

#### Entrée

* Ensemble des fichiers TSV produits par `extract_xy_intensities`.

#### Sortie intermédiaire

Fichier :

* `all_xy_intensities.tsv`

contenant les intensités des chromosomes X et Y pour l’ensemble des échantillons présents dans la `sample_sheet`.

Ce fichier n’est pas exporté dans le répertoire final des résultats. Il constitue un fichier intermédiaire utilisé par les étapes suivantes du workflow.

### `control_sex`

Ce process réalise un contrôle d’identitovigilance basé sur le sexe des échantillons en comparant le sexe prédit à partir des intensités des chromosomes X et Y avec l’information renseignée dans la `sample_sheet`.

#### Entrée

* **xy_intensities_merged** : fichier produit par le process `merge_xy_intensities`.

#### Sortie

##### `control_sex_report`

Fichier PDF contenant :

* un graphique représentant chaque échantillon selon ses intensités sur les chromosomes X et Y ;
* l’indication de la concordance entre le sexe prédit et le sexe renseigné dans la `sample_sheet` ;
* un tableau récapitulatif indiquant, pour chaque échantillon affiché, son statut :

  * `correct`
  * `mismatch`
  * `undetermined`

Lorsque le nombre d’échantillons est supérieur à 24, seuls les échantillons en échec sont affichés dans le tableau.

##### `all_predicted_sex.tsv`

Fichier TSV contenant, pour l’ensemble des échantillons de la `sample_sheet` :

* les intensités des chromosomes X et Y ;
* le sexe renseigné ;
* le sexe prédit ;
* le statut de concordance entre ces informations.


## Module CNV 

Le module CNV permet d’inférer les variations du nombre de copies (CNV) globales et focales à partir des données de méthylation issues de puces Illumina.

L’inférence des CNV repose sur le package R *conumee2*, qui constitue actuellement l’une des méthodes de référence pour l’analyse des CNV à partir de données de méthylation. Cette méthode est décrite dans la publication :

> Daenekas B. et al., *Conumee 2.0: enhanced copy-number variation analysis from DNA methylation arrays for humans and mice*, Bioinformatics, 2024.

Le package *conumee2* présente notamment les caractéristiques suivantes :

* compatibilité avec les puces 450k, EPICv1 et EPICv2 ;
* intégration avec les principaux objets Bioconductor ;
* normalisation adaptée aux données de méthylation (*tangent normalization*) ;
* détection optionnelle d’altérations focales par bootstrap ;
* génération automatisée de profils et de visualisations CNV.

Les CNV inférés à partir des données de méthylation constituent une estimation indirecte des altérations génomiques. Leur résolution dépend notamment de la densité et de la répartition des sondes, du bruit expérimental, du jeu de données de référence utilisé ainsi que des paramètres de binning et de segmentation.

Le workflow mis en œuvre par *conumee2* est le suivant :

1. combinaison des intensités méthylées et non méthylées ;
2. normalisation par rapport à un jeu d’échantillons de référence (*tangent normalization*) ;
3. regroupement des signaux en bins génomiques ;
4. segmentation des profils afin d’identifier les régions présentant des gains ou pertes de copies ;
5. détection optionnelle d’altérations focales par bootstrap ;
6. génération des profils et visualisations CNV.

![Workflow du package conumee2](https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/bioinformatics/40/2/10.1093_bioinformatics_btae029/2/btae029f1.jpeg?Expires=1787322780&Signature=jMKkjU86p9pmGwRLgd3l4BXIwr83bJJL3Hfj1XR-I9N8rKxtrt1S2R2Nj7hAGYBX0gY91UihHhtY4jS01dzMRwEST3ca6Z~ZSPHDvqjEZC7V1n3nzVeoePg6PrmWpWXmggIEX84Jx02abTmDMGLtxuJGxU8bluUOZ0~CgYxx08imcSDiBiY4hevb5FUvooaQIPU95hWKbOWL1TM~UaOrpTfGMLHQx0V07mKCFpMH25~S2XVezOlWVjEQiRVPsYfK1MDmk8RVwMSG1LHOd9Uhc6gtcPU5WmBN79If~PvDk7pZFtSMDV~s51wupTgdZJ7cjxi-6TcWLAN9jqpVqi2Dgw__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA)

*Workflow du package R conumee2*
### Gestion des chromosomes sexuels

Les chromosomes sexuels présentent des profils de méthylation et de dosage susceptibles d’influencer la normalisation et la segmentation.

Afin de limiter ces effets, deux jeux de données de référence distincts sont utilisés :

* une référence masculine ;
* une référence féminine ; 
* une référence mixte.

Chaque échantillon est analysé avec la référence correspondant au sexe renseigné dans la `sample_sheet`. Lorsque le sexe renseigné est "U", la référence mixte est utilisée.

### Particularité des puces EPICv2

Lors de l’utilisation de *conumee2* avec des données issues de puces *Infinium MethylationEPIC v2.0*, la fonction `CNV.fit` ne fonctionne pas lorsqu’un seul échantillon est chargé.

Afin de contourner cette limitation, l'échantillon est duppliqué puis renommé afin de passer en entrée de `CNV.fit` deux échantillons en EPICv2.

### Jeu de données de référence

Le choix du jeu de données de référence influence directement :

* le niveau de bruit résiduel après normalisation ;
* la stabilité de la segmentation ;
* la détection des altérations focales ;
* les seuils de significativité des segments.

Les références fournies par défaut dans le pipeline sont construites à partir d’échantillons non tumoraux du dépôt GEO [GSE306226](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE306226). Seuls les neurones et les microglies ont été conservés afin de constituer un jeu de référence composé de tissus sains.

### Recommandations

#### Construction d’un jeu de référence

Pour construire un jeu de données de référence adapté, les recommandations suivantes sont appliquées :

* utiliser au minimum 16 échantillons contrôles ;
* utiliser le même type tissulaire ;
* utiliser le même protocole expérimental ;
* inclure des échantillons FFPE et frozen lorsque les deux types doivent être analysés.

#### Qualité des données

Le bruit constitue un indicateur important de la qualité des profils CNV.

Il est calculé selon la formule suivante :

```text
noise = sqrt( (1/(n-1)) * sum(R_{i+1} - R_i)^2 )
```

Lorsque cette valeur est supérieure à `0.6`, l’exclusion de l’échantillon est recommandée.

#### Paramètres de binning

Les paramètres de binning utilisés par défaut correspondent aux paramètres initialement optimisés pour les puces *Infinium Human Methylation 450K BeadChip*.

Ils représentent un compromis entre résolution et niveau de bruit. La taille des bins est fixée par défaut à environ 50 kb. Ces paramètres peuvent être adaptés selon le type de puce utilisé et l’objectif de l’analyse.

### (`compute_cnv`)

Ce process calcule les CNV pour chaque échantillon de la `sample_sheet` et fournit des visualisations et des fichiers tabulés portant ces informations.

#### Entrée

Le process `compute_cnv` nécessite les éléments suivants :

* Tuple issu du process `load_idats_minfi` :

  * **rgset_rds** ; 
  * **sex** ;
  * **sample_id**.
* **ref_f** : objet `.Rdata` contenant la référence féminine ;
* **ref_m** : objet `.Rdata` contenant la référence masculine ;
* **ref_nm** : objet `.Rdata` contenant la référence mixte ;
* **anno** : objet `.Rdata` contenant l’annotation de genomic binning ;
* **CNV_focal** : True ou False en groovy déterminant le calcul ou non des CNV focaux et donc de la dépendance à une connexion internet.

#### Sortie

Pour chaque échantillon, les résultats sont organisés sous forme de figures et de fichiers de données.

##### Visualisations

* `sample_id_AllChr.png` : profil CNV global comprenant les chromosomes 1 à 22 ainsi que les chromosomes sexuels ;
* `sample_id_chrXX.png` : profils CNV individuels par chromosome ;
* `sample_id_detailed_region.png` : visualisations des régions d’intérêt ;
* `sample_id_Genes.png` : synthèse des CNV pour l’ensemble des régions d’intérêt.

#### Fichiers de données

* `sample_id_CNVsegments.tsv` : fichier TSV contenant les segments CNV et leurs p-values associées ;
* `sample_id_CNVbins.igv` : fichier IGV contenant les valeurs agrégées par bin ;
* `sample_id_CNVdetail.tsv` : fichier TSV contenant les valeurs des régions d’intérêt ;
* `sample_id_CNVprobes.igv` : fichier IGV contenant les valeurs au niveau des sondes ;
* `sample_id_metrics.tsv` : fichier TSV contenant les valeurs de bruit, de shift et les coefficients de la régression liée au calcul des CNV.

Les résultats sont répartis dans différents répertoires selon leur nature :

* `img` pour les visualisations ;
* `tables` pour les fichiers tabulaires.

Les visualisations par chromosome et par région d’intérêt sont regroupées dans des sous-répertoires dédiés.


### `merge_cnv_detail` 

Ce process concatène les fichiers **detail** générés par `compute_cnv` en un unique fichier TSV.

#### Entrée 

Le process `merge_cnv_detail`nécessite : 

* all_CNV_detail ;
* l'ensemble des fichiers detail.TXT produits par `compute_cnv`.

#### Sortie

Fichier :

* `all_CNV_detail.tsv`

contenant les valeurs des régions d’intérêt (gènes) pour l’ensemble des échantillons présents dans la `sample_sheet`.


### `merge_cnv_segment` 

Ce process concatène les fichiers **segment** générés par `compute_cnv` en un unique fichier TSV.

#### Entrée 

Le process `merge_cnv_segment`nécessite : 

* all_CNV_segment ;
* l'ensemble des fichiers SEG produits par `compute_cnv`.

#### Sortie

Fichier :

* `all_CNV_segment.tsv`

contenant les segments CNV et leurs p-values associées pour l’ensemble des échantillons présents dans la `sample_sheet`.


### `merge_cnv_metrics` 

Ce process concatène les fichiers **metrics** générés par `compute_cnv` en un unique fichier TSV.

#### Entrée 

Le process `merge_cnv_metrics`nécessite : 

* all_CNV_metrics ;
* l'ensemble des fichiers metrics.txt produits par `compute_cnv`.

#### Sortie

Fichier :

* `all_CNV_metrics.tsv`

contenant les valeurs de bruit, de shift et les coefficients de la régression liée au calcul des CNV pour l’ensemble des échantillons présents dans la `sample_sheet`.


## Module MGMT

### Description

Le module MGMT permet de prédire l’état de méthylation du promoteur du gène *MGMT* à partir des profils de méthylation issus des puces Illumina.

La prédiction repose sur le package R `mgmtstp27`, qui utilise les niveaux de méthylation des sondes **cg12434587** et **cg12981137** afin d’estimer le statut MGMT de l’échantillon.

Le statut de méthylation de *MGMT* constitue un biomarqueur moléculaire largement utilisé dans les gliomes. La méthylation du promoteur de ce gène est associée à une diminution de son expression et à une meilleure sensibilité aux agents alkylants tels que le témozolomide.

Dans ce module, aucune normalisation n’est appliquée avant la prédiction réalisée par le package `mgmtstp27`. Les calculs sont effectués directement à partir des données chargées par le workflow.

### `mgmt`

#### Entrée

Tuple issu du process `load_idats_minfi` :

* **rgset_rds** ; 
* **sex** ;
* **sample_id**.

#### Sortie

Pour chaque échantillon :

* fichier TSV `sample_id_mgmt.tsv` contenant les valeurs e le statut de méthylation du promoteur du gène *MGMT* de l'échantillon ;

* fichier PDF `MGMT_plot_minfi_sample_id.pdf`.

Ce rapport contient :

* une représentation graphique de la prédiction du statut *MGMT* ;
* le statut prédit (*methylated* ou *unmethylated*) ;
* un tableau récapitulatif des intervalles de confiance associés à la prédiction, calculés sous l’hypothèse d’une distribution normale.

### `merge_mgmt` 

Ce process concatène les fichiers TSV générés par `mgmt` en un unique fichier TSV.

#### Entrée 

Le process `merge_mgmt`nécessite : 

* all_mgmt ;
* l'ensemble des fichiers mgmt.tsv produits par `mgmt`.

#### Sortie

Fichier :

* `all_mgmt.tsv`

contenant les valeurs et le statut de méthylation du promoteur du gène *MGMT* pour l’ensemble des échantillons présents dans la `sample_sheet`.


## Module Tumor Purity

### Description

Le module de pureté tumorale permet d’estimer la proportion de cellules tumorales présentes dans un échantillon à partir des données de méthylation.

L’estimation repose sur le package R `RFpurify`, qui implémente deux modèles de régression par *Random Forest* entraînés à reproduire deux méthodes de référence de mesure de la pureté tumorale :

* **ABSOLUTE**, basée sur les altérations du nombre de copies ;
* **ESTIMATE**, basée sur les signatures transcriptomiques stromales et immunitaires.

Contrairement à certaines approches nécessitant des échantillons contrôles ou une connaissance préalable du type tumoral, les modèles utilisés par `RFpurify` permettent d’estimer la pureté tumorale directement à partir des profils de méthylation.

Les modèles fournis par `RFpurify` ont été entraînés à partir de données issues de puces **Infinium Human Methylation 450K BeadChip**. Lors de l’analyse de données provenant de puces **Infinium MethylationEPIC v2.0**, plusieurs adaptations sont nécessaires :

* certaines sondes présentes sur les puces 450k ont été renommées sur les puces EPICv2 ;
* certaines sondes utilisées par les modèles de prédiction ne sont plus présentes sur les puces EPICv2.

Afin de conserver la compatibilité avec les modèles de `RFpurify`, les nouveaux identifiants de sondes sont convertis vers leurs identifiants historiques.

Pour les sondes absentes des puces EPICv2 mais requises par les modèles, une valeur arbitraire de 0.5 est attribuée. Cette procédure a été testée et validée afin de prendre en compte l’incertitude induite par l’absence de ces sondes.

Les valeurs rapportées correspondent alors :

* à la prédiction obtenue par la méthode ABOSULTE ;
* à la prédiction obtenue par la méthode ESTIMATE.

### `compute_tumor_purity`

#### Entrée

Tuple issu du process `load_idats_minfi` :

* **rgset_rds** ;
* **sex** ;
* **sample_id**.

#### Sortie

Pour chaque échantillon :

* fichier TSV `sample_id_tumor_purity.tsv` contenant les valeurs de pureté tumorale prédites par les méthodes ABSOLUTE et ESTIMATE ; 

* fichier PDF `sample_id_tumor_purity.pdf`.

Ce rapport contient :

* une représentation graphique des estimations de pureté tumorale obtenues par les méthodes `absolute` et `estimate` au cours des 100 itérations ;
* un tableau récapitulatif présentant :

  * la moyenne des valeurs prédites ;
  * l’intervalle de confiance associé à chaque méthode.

### `merge_tumor_purity` 

Ce process concatène les fichiers TSV générés par `tumor_purity` en un unique fichier TSV.

#### Entrée 

Le process `merge_tumor_purity`nécessite : 
* all_tumor_purity ;
* l'ensemble des fichiers tumor_purity.tsv produits par `tumor_purity`.

#### Sortie

Fichier :

* `all_tumor_purity.tsv`

contenant les valeurs prédites par les méthodes ABSOLUTE et ESTIMATE pour l’ensemble des échantillons présents dans la `sample_sheet`.


## Ressources supplémentaires 

### Où trouver les informations d'annotation dans les différents objets : 

1. **ewastools**

- *Sample_IDAT* :

`x[["meta"]][["sample_id"]]`

- *Probe_ID* :

`x[["manifest"]][["ilmn_id"]]`

- *Probe_design* :

`x[["manifest"]][["probe_design"]]`

- *Nom des sondes contrôle* :
 
`x[["controls"]][["group"]]`

- *Valeurs des sondes de contrôle* :

`x[["ctrlR"]]`
`x[["ctrlG"]]`

- *Platform* :

`x[["platform"]]`

2. **minfi : rgset**

- *Sample_IDAT* :

`x@colData@rownames`

- *Probe_ID* :

```R
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
rownames(getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38))
```

- *Annotation* :

`getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)$Type`

3. **conumee2**

- *Sample_IDAT* : 

`names(x[i])`


### Preprocessing

#### Normalisation avec `minfi`

Les méthodes de normalisation du package R `minfi` s’appliquent sur des objets de type `RGset`, qui sont ensuite convertis en objets `Mset`.

Plusieurs méthodes de preprocessing sont implémentées :

* `preprocessRaw`  
  → aucune normalisation appliquée ;

* `preprocessIllumina`  
  → méthode de normalisation utilisée par Illumina GenomeStudio ;

* `preprocessNoob`  
  → méthode *normal-exponential out-of-band* permettant une correction du bruit de fond et du biais de coloration ;

* `preprocessQuantile`  
  → méthode de *stratified quantile normalization* pour les puces de méthylation Illumina. Les sondes sont stratifiées selon leur région génomique (CpG islands, shores, etc.).

Une méthode de normalisation modifiée est également disponible dans le dépôt GitHub associé à la publication du classifier v12.8 du DKFZ : `MNPpreprocessIllumina`.

Cette méthode repose sur `preprocessIllumina` mais réalise une normalisation à 10 000 au lieu d’utiliser l’intensité moyenne des sondes de contrôle de normalisation du premier échantillon/de l’échantillon de référence.


#### Normalisation avec `ewastools`

Les méthodes de normalisation du package R `ewastools` s’appliquent sur les objets issus de la fonction `read_idats`.

Plusieurs types de traitements sont possibles selon le type d’objet de sortie souhaité.

##### **Objets bruts** (*raw objects*)

* `correct_dye_bias`
* `correct_dye_bias2`

##### **Beta-values**

* `dont_normalize`
* `normalize`

Ces deux types de normalisation peuvent être combinés : 

```R
beta = meth %>% correct_dye_bias %>% dont_normalize
```
