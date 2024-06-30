# Changelog

## [0.2.0](https://github.com/bihealth/snappy-pipeline/compare/v0.1.0...v0.2.0) (2024-06-28)


### Features

* use pydantic models for validation ([#496](https://github.com/bihealth/snappy-pipeline/issues/496)) ([bf39678](https://github.com/bihealth/snappy-pipeline/commit/bf3967854337522484417d8fc6e1e0ac67f0c43d))

## 0.1.0 (2024-06-28)


### ⚠ BREAKING CHANGES

* pin dependency version in all wrapper conda env yamls ([#492](https://github.com/bihealth/snappy-pipeline/issues/492))
* broken plot generation for the control_freec tool in the somatic_wgs_cnv_calling step, and the sample selection for the mutect2 panel_of_normal generation can be based on a user-generated file.
* use snakemake tmpdir resource ([#319](https://github.com/bihealth/snappy-pipeline/issues/319)) (#320)
* merge wgs_{sv,cnv,mei}_calling into sv_calling_wgs ([#275](https://github.com/bihealth/snappy-pipeline/issues/275)) (#315)
* merge {variant,wgs_{sv,cnv}}_calling_export into varfish_export ([#308](https://github.com/bihealth/snappy-pipeline/issues/308)) (#309)
* move "bcftools roh" as report in variant_calling ([#310](https://github.com/bihealth/snappy-pipeline/issues/310)) (#311)
* rename targeted_seq_cnv_calling to sv_calling_targeted ([#305](https://github.com/bihealth/snappy-pipeline/issues/305)) (#306)
* rename helper_gcnv_model_target_seq to helper_gcnv_targeted ([#303](https://github.com/bihealth/snappy-pipeline/issues/303)) (#304)

### Features

* (430) picard metrics for bam files ([#431](https://github.com/bihealth/snappy-pipeline/issues/431)) ([35b24cd](https://github.com/bihealth/snappy-pipeline/commit/35b24cd70e4352244b6d09d414c6fae761c7dece))
* 237 revive cnvkit ([#339](https://github.com/bihealth/snappy-pipeline/issues/339)) ([bc4c55b](https://github.com/bihealth/snappy-pipeline/commit/bc4c55b2b095a9a8c792cf937efd6fda87d8e09e))
* 385 cbioportal export wrapper missing information ([#411](https://github.com/bihealth/snappy-pipeline/issues/411)) ([0ac31e8](https://github.com/bihealth/snappy-pipeline/commit/0ac31e861af0edfe67c1b82d086298b2d9fc9081))
* 438 vep annotation order ([#441](https://github.com/bihealth/snappy-pipeline/issues/441)) ([c28d67e](https://github.com/bihealth/snappy-pipeline/commit/c28d67e23246d319fa88a20331d2324eee379a9f))
* add CI for a scaled down version of a cancer WES pipeline ([#499](https://github.com/bihealth/snappy-pipeline/issues/499)) ([b971088](https://github.com/bihealth/snappy-pipeline/commit/b9710883206aa9629c357940acf24bd03493fcd3))
* add delly2 and manta for sv_calling_targeted ([#277](https://github.com/bihealth/snappy-pipeline/issues/277)) ([#307](https://github.com/bihealth/snappy-pipeline/issues/307)) ([99dae5f](https://github.com/bihealth/snappy-pipeline/commit/99dae5ffc1b828b7e7b46e1fa2259223e3a5c243))
* add missense TMB calculation ([#432](https://github.com/bihealth/snappy-pipeline/issues/432)) ([4fd0c73](https://github.com/bihealth/snappy-pipeline/commit/4fd0c73c4bd75bf4aa5ef9a464713a2c5d4fada1))
* Added support for scarHRD ([#429](https://github.com/bihealth/snappy-pipeline/issues/429)) ([c5f9dde](https://github.com/bihealth/snappy-pipeline/commit/c5f9dde2d2ce8570779030c1917b3132918ad33b))
* added vep to the somatic variant annotation tools ([#384](https://github.com/bihealth/snappy-pipeline/issues/384)) ([35cbe09](https://github.com/bihealth/snappy-pipeline/commit/35cbe0974ad4d56419ae8790110e418c22ef84f7))
* adding bwa_mem2 wrapper ([#279](https://github.com/bihealth/snappy-pipeline/issues/279)) ([#312](https://github.com/bihealth/snappy-pipeline/issues/312)) ([1374418](https://github.com/bihealth/snappy-pipeline/commit/13744188bb1943b007ed3ad02a6ab6367c865fee))
* Adding tmb calculation  ([#403](https://github.com/bihealth/snappy-pipeline/issues/403)) ([706a9f5](https://github.com/bihealth/snappy-pipeline/commit/706a9f561fe96a9469fec15f232940b265266e8e))
* allow helper_gcnv_model_target_seq to build for multiple kits at the same time ([#285](https://github.com/bihealth/snappy-pipeline/issues/285)) ([#293](https://github.com/bihealth/snappy-pipeline/issues/293)) ([efe6539](https://github.com/bihealth/snappy-pipeline/commit/efe653906ec841496f9b1a8c32cfa7162b7fe242))
* also write out mapping quality bigWig file ([#476](https://github.com/bihealth/snappy-pipeline/issues/476)) ([#480](https://github.com/bihealth/snappy-pipeline/issues/480)) ([dc1d92b](https://github.com/bihealth/snappy-pipeline/commit/dc1d92b729850ac6571b1f3f9a43b9d61ff531a2))
* annotating SV coverage with varfish-annotator ([#371](https://github.com/bihealth/snappy-pipeline/issues/371)) ([e9fb880](https://github.com/bihealth/snappy-pipeline/commit/e9fb8807a1c48854ce28e8b7942e0f3b4d4a7765))
* bringing back variant_annotation with VEP ([#386](https://github.com/bihealth/snappy-pipeline/issues/386)) ([#387](https://github.com/bihealth/snappy-pipeline/issues/387)) ([d4ac11e](https://github.com/bihealth/snappy-pipeline/commit/d4ac11ea5ad84fad6b8c9ff5db69d0dd446b8fef))
* bump and fix Mehari for SVs in varfish_export ([#444](https://github.com/bihealth/snappy-pipeline/issues/444)) ([0ca23f6](https://github.com/bihealth/snappy-pipeline/commit/0ca23f645d38130ac08f85b1861ae795f62f96b2))
* bump ngs-chew to 0.6.x, adding aafs to output ([#317](https://github.com/bihealth/snappy-pipeline/issues/317)) ([feb13a3](https://github.com/bihealth/snappy-pipeline/commit/feb13a352116ce5fbe7123f9aafbe1f1122026f1))
* bump ngs-chew to 0.7 ([#334](https://github.com/bihealth/snappy-pipeline/issues/334)) ([7f1fa1a](https://github.com/bihealth/snappy-pipeline/commit/7f1fa1ae3dec9416df9c8e35e542f67df6e09811))
* bump ngs-chew to v0.9.x ([#478](https://github.com/bihealth/snappy-pipeline/issues/478)) ([c54bba1](https://github.com/bihealth/snappy-pipeline/commit/c54bba11a558ad1cb7dfa5f785b4f7b7aa526ea2))
* bump varfish-annotator to v0.32 ([#365](https://github.com/bihealth/snappy-pipeline/issues/365)) ([3239564](https://github.com/bihealth/snappy-pipeline/commit/32395640f320e3815fa0bc7953fb3a2eecff99a2))
* bump varfish-annotator-cli to v0.30 ([#314](https://github.com/bihealth/snappy-pipeline/issues/314)) ([#336](https://github.com/bihealth/snappy-pipeline/issues/336)) ([3d045e6](https://github.com/bihealth/snappy-pipeline/commit/3d045e6ad791823b3e6b28f5cf3e2b541c70bf17))
* change varfish-annotator to mehari ([#392](https://github.com/bihealth/snappy-pipeline/issues/392)) ([#393](https://github.com/bihealth/snappy-pipeline/issues/393)) ([c86c665](https://github.com/bihealth/snappy-pipeline/commit/c86c665b8132e462c7916c67c171192abaff4a43))
* feature-complete end-to-end run of the Becnel dataset on hs37d5 & GRCh38.d1.vd1 ([#460](https://github.com/bihealth/snappy-pipeline/issues/460)) ([4874074](https://github.com/bihealth/snappy-pipeline/commit/4874074d8200ac33d29377ce1a812dfbb48adc31))
* Fixing Mantis ([#427](https://github.com/bihealth/snappy-pipeline/issues/427)) ([e008e7e](https://github.com/bihealth/snappy-pipeline/commit/e008e7edc42371eaee00e3a3861816f2ddfdb39c))
* implement ngs_mapping fingerprinting ([#283](https://github.com/bihealth/snappy-pipeline/issues/283)) ([#289](https://github.com/bihealth/snappy-pipeline/issues/289)) ([a3303ac](https://github.com/bihealth/snappy-pipeline/commit/a3303ac4964c988c88c5342712046189cab77183))
* Initial implementation of Molecular bar codes handling using AGeNT ([#462](https://github.com/bihealth/snappy-pipeline/issues/462)) ([768dded](https://github.com/bihealth/snappy-pipeline/commit/768dded3bdfc3d86a1029b88bf2e25c224eee8ee))
* integrate GATK4 HaplotypeCaller ([#287](https://github.com/bihealth/snappy-pipeline/issues/287)) ([#299](https://github.com/bihealth/snappy-pipeline/issues/299)) ([2aeaa01](https://github.com/bihealth/snappy-pipeline/commit/2aeaa01249dd79e654acb5455bf49221f0c211a5))
* make germline variants pon optional for mutect2 ([#375](https://github.com/bihealth/snappy-pipeline/issues/375)) ([30bc591](https://github.com/bihealth/snappy-pipeline/commit/30bc591cdf23ebc5ffe288e6bcd76e00311256d3))
* make ngs_mapping sub steps generate links ([#291](https://github.com/bihealth/snappy-pipeline/issues/291)) ([#294](https://github.com/bihealth/snappy-pipeline/issues/294)) ([6f2edd7](https://github.com/bihealth/snappy-pipeline/commit/6f2edd74f423ed7c65776e4f3943f228fca79990))
* melt for sv_calling_targeted ([#321](https://github.com/bihealth/snappy-pipeline/issues/321)) ([13ee84b](https://github.com/bihealth/snappy-pipeline/commit/13ee84b60bdf3e168d8d3dea55425314decf59f6))
* merge {variant,wgs_{sv,cnv}}_calling_export into varfish_export ([#308](https://github.com/bihealth/snappy-pipeline/issues/308)) ([#309](https://github.com/bihealth/snappy-pipeline/issues/309)) ([7c35ec3](https://github.com/bihealth/snappy-pipeline/commit/7c35ec38d6f861b5ee53e5afed8ccffc3c78dce0))
* merge gCNV output in case of multi-kit families ([#292](https://github.com/bihealth/snappy-pipeline/issues/292)) ([#465](https://github.com/bihealth/snappy-pipeline/issues/465)) ([ee64a98](https://github.com/bihealth/snappy-pipeline/commit/ee64a98a6c29d307b997209d66b5289e39e88621))
* ngs_mapping cleanup, versions bump ([#278](https://github.com/bihealth/snappy-pipeline/issues/278)) ([#280](https://github.com/bihealth/snappy-pipeline/issues/280)) ([af8e21f](https://github.com/bihealth/snappy-pipeline/commit/af8e21f264ab93f143d4e7d5b89ee1d99889f247))
* parallel execution in gatk4_hc_gvcf combine_gvcfs and genotype ([#318](https://github.com/bihealth/snappy-pipeline/issues/318)) ([4e888e0](https://github.com/bihealth/snappy-pipeline/commit/4e888e068e96d81c646112fb5d16305fa0ea5ab8))
* path_exon_bed can be provided independent of data type ([#389](https://github.com/bihealth/snappy-pipeline/issues/389)) ([#390](https://github.com/bihealth/snappy-pipeline/issues/390)) ([7759d4d](https://github.com/bihealth/snappy-pipeline/commit/7759d4dfa9487ecae72fc4841b8b086ecb7bec0d))
* PureCN implementation ([#453](https://github.com/bihealth/snappy-pipeline/issues/453)) ([e394e8c](https://github.com/bihealth/snappy-pipeline/commit/e394e8cd4bd76eef4725bbaa9ff7711ab5aa9dcc))
* remove legacy variant callers ([#295](https://github.com/bihealth/snappy-pipeline/issues/295)) ([#297](https://github.com/bihealth/snappy-pipeline/issues/297)) ([fa2fc6c](https://github.com/bihealth/snappy-pipeline/commit/fa2fc6c4e8fc2a9797325a09fa06401d43e35cd9))
* remove plot-bamstats ([#288](https://github.com/bihealth/snappy-pipeline/issues/288)) ([#290](https://github.com/bihealth/snappy-pipeline/issues/290)) ([125b601](https://github.com/bihealth/snappy-pipeline/commit/125b6013892f502bbce336403a821395e60b16b0))
* replacing custom coverage QC scripts with "alfred qc" ([#481](https://github.com/bihealth/snappy-pipeline/issues/481)) ([#482](https://github.com/bihealth/snappy-pipeline/issues/482)) ([ae8fe98](https://github.com/bihealth/snappy-pipeline/commit/ae8fe9859a54690ed4908f8d1a20dd57418f96e0))
* Somatic cnv checking ([#426](https://github.com/bihealth/snappy-pipeline/issues/426)) ([18217bc](https://github.com/bihealth/snappy-pipeline/commit/18217bc72a36b2a7dc35cecf2cc3a2d07456f836))
* support for snakemake --batch, pipeline.sh for array jobs ([#325](https://github.com/bihealth/snappy-pipeline/issues/325)) ([#326](https://github.com/bihealth/snappy-pipeline/issues/326)) ([6ffa5cc](https://github.com/bihealth/snappy-pipeline/commit/6ffa5ccb8edad9f4d73aada1feca01c933844f8e))
* support for somatic variant calling without normal ([#503](https://github.com/bihealth/snappy-pipeline/issues/503)) ([66f555c](https://github.com/bihealth/snappy-pipeline/commit/66f555c40caf5b4b1e2e4c17bf6ffc8169ba8baf))
* switch GATK 3 HC/UG to GNU parallel based multithreading ([#301](https://github.com/bihealth/snappy-pipeline/issues/301)) ([#302](https://github.com/bihealth/snappy-pipeline/issues/302)) ([e23e300](https://github.com/bihealth/snappy-pipeline/commit/e23e3006da2b3167746162f3ef072ccc4c6d9b45))
* upgrade to GATK v4.3.0 for gCNV variant calling ([#218](https://github.com/bihealth/snappy-pipeline/issues/218)) ([#265](https://github.com/bihealth/snappy-pipeline/issues/265)) ([372464e](https://github.com/bihealth/snappy-pipeline/commit/372464eac8bdb9500499676fdb735bc75ad54e26))
* use filtered variants in tmb and signatures ([#470](https://github.com/bihealth/snappy-pipeline/issues/470)) ([3e3cead](https://github.com/bihealth/snappy-pipeline/commit/3e3cead3a9ea02484a0d7d065722282ee3f4f112))
* use snakemake tmpdir resource ([#319](https://github.com/bihealth/snappy-pipeline/issues/319)) ([#320](https://github.com/bihealth/snappy-pipeline/issues/320)) ([c8512b0](https://github.com/bihealth/snappy-pipeline/commit/c8512b0b9c7f3eaece974973a916b96eee65dd01))


### Bug Fixes

* 355 wgs svcnv export external are also affected by changes in the wrapper for varfish annotator ([#356](https://github.com/bihealth/snappy-pipeline/issues/356)) ([a660bda](https://github.com/bihealth/snappy-pipeline/commit/a660bda9ca3a7478b68f9ad6c16c9f7a49f18cc1))
* Add bwa to mbcs meta-tool environment, closes [#518](https://github.com/bihealth/snappy-pipeline/issues/518) ([b3305c5](https://github.com/bihealth/snappy-pipeline/commit/b3305c5bf918b10b9b1240a4108e08d2c03383f8))
* add missing mask_duplicates to bwa_mem2 wrapper ([#324](https://github.com/bihealth/snappy-pipeline/issues/324)) ([c5cc7a0](https://github.com/bihealth/snappy-pipeline/commit/c5cc7a0504f0ed06565331cfe61f343ec46a42de))
* adjust STAR to new ngs_mapping ([#366](https://github.com/bihealth/snappy-pipeline/issues/366)) ([#367](https://github.com/bihealth/snappy-pipeline/issues/367)) ([65c5b6e](https://github.com/bihealth/snappy-pipeline/commit/65c5b6ee683a84c06e502af4ed185e4927118589))
* adjustment of manta wrapper for sv_calling_wgs ([#351](https://github.com/bihealth/snappy-pipeline/issues/351)) ([e2cd197](https://github.com/bihealth/snappy-pipeline/commit/e2cd19700e55e1bfc3936a0d5b9a16ed5ddcdc4f))
* allow to pass single batch number to pipeline_job.sh ([#327](https://github.com/bihealth/snappy-pipeline/issues/327)) ([#333](https://github.com/bihealth/snappy-pipeline/issues/333)) ([1c09320](https://github.com/bihealth/snappy-pipeline/commit/1c09320373b8ab84d65f0f3fe611844fea5705e4))
* Always include alfred qc files in input for varfish export bam QC [#511](https://github.com/bihealth/snappy-pipeline/issues/511) ([#516](https://github.com/bihealth/snappy-pipeline/issues/516)) ([d13e595](https://github.com/bihealth/snappy-pipeline/commit/d13e59526dd382c51a3065d0c1e6ca68b59a92f7))
* bump mehari to 0.21.1 fixes effect prediction ([#474](https://github.com/bihealth/snappy-pipeline/issues/474)) ([c0ba296](https://github.com/bihealth/snappy-pipeline/commit/c0ba2964b4c28bbef154f65fcb7afed324ffe3ca))
* bump mehari to v0.25.5 for bugs fixed there ([#510](https://github.com/bihealth/snappy-pipeline/issues/510)) ([95b451b](https://github.com/bihealth/snappy-pipeline/commit/95b451b1b259b24c29b0f9453a521fcc2f1b92ba))
* bump melt memory further ([#374](https://github.com/bihealth/snappy-pipeline/issues/374)) ([fad1ad6](https://github.com/bihealth/snappy-pipeline/commit/fad1ad6ec184882ef138ee186a2f9a0cf6a11cf5))
* bump time of ngs_mapping/bam_collect_doc ([#394](https://github.com/bihealth/snappy-pipeline/issues/394)) ([04551d5](https://github.com/bihealth/snappy-pipeline/commit/04551d5ed26b31e6c822df1a1fe6d83587edec74))
* bwa_mem2 wrapper only used right read ([#316](https://github.com/bihealth/snappy-pipeline/issues/316)) ([53c88aa](https://github.com/bihealth/snappy-pipeline/commit/53c88aaf784f07125e697d35faf330f58186613e))
* check whether gCNV is enabled correctly ([#357](https://github.com/bihealth/snappy-pipeline/issues/357)) ([e5ec1b3](https://github.com/bihealth/snappy-pipeline/commit/e5ec1b367f6bcf0d86be45de093e488f9b255ac2))
* computation of MD5 sums in gatk4_hc genotype wrapper ([#338](https://github.com/bihealth/snappy-pipeline/issues/338)) ([0900d9d](https://github.com/bihealth/snappy-pipeline/commit/0900d9d5c4fa96230a7c622d118e86ed19b8e5a8))
* conda env for merge_multikit_families - install attrs package over attr ([#487](https://github.com/bihealth/snappy-pipeline/issues/487)) ([fab3e36](https://github.com/bihealth/snappy-pipeline/commit/fab3e3687e165acf76114b1a025a07c2542c3d2f))
* correct name of annotate_svs in Snakefiles ([#344](https://github.com/bihealth/snappy-pipeline/issues/344)) ([#345](https://github.com/bihealth/snappy-pipeline/issues/345)) ([f1355ae](https://github.com/bihealth/snappy-pipeline/commit/f1355ae88673bc291b938e7806b3838bf5abdbe5))
* definition of gcnv conda enviroment for mamba&gt;=1.3 ([#452](https://github.com/bihealth/snappy-pipeline/issues/452)) ([#454](https://github.com/bihealth/snappy-pipeline/issues/454)) ([50d0781](https://github.com/bihealth/snappy-pipeline/commit/50d0781373014eec97b2e993bb630ed922bb0847))
* do not expect md5 sums for make_vcf step of melt ([#323](https://github.com/bihealth/snappy-pipeline/issues/323)) ([65b45bd](https://github.com/bihealth/snappy-pipeline/commit/65b45bdc5c9dd197c32c5ed98668c639b4c3e0b8))
* eb_filter wrapper assuming ANN/CSQ is present ([#506](https://github.com/bihealth/snappy-pipeline/issues/506)) ([4dc4948](https://github.com/bihealth/snappy-pipeline/commit/4dc49488cc0d04899e2202a6c6613d852c1d3197))
* enforce ploidy in GATK gCNV ploidy calling wrapper ([#377](https://github.com/bihealth/snappy-pipeline/issues/377)) ([f244691](https://github.com/bihealth/snappy-pipeline/commit/f2446915e4f4b7da25ac09d74005b366bc1806f3))
* explicitely enable conda when mamba is disabled ([#343](https://github.com/bihealth/snappy-pipeline/issues/343)) ([20bb62c](https://github.com/bihealth/snappy-pipeline/commit/20bb62cd2ae7fc7ed8f795133946b7a709504363))
* explicitly require gawk in mutect2/filter wrapper environment ([#498](https://github.com/bihealth/snappy-pipeline/issues/498)) ([6993233](https://github.com/bihealth/snappy-pipeline/commit/6993233bc8b50641860d8a4bfc70b21a4ff09076))
* fix somatic snappy ([#477](https://github.com/bihealth/snappy-pipeline/issues/477)) ([#491](https://github.com/bihealth/snappy-pipeline/issues/491)) ([af961cc](https://github.com/bihealth/snappy-pipeline/commit/af961cca3ef4d8c7a5a4a2ec08c3161260ce9f5e))
* fix trimadap-mt & bbduk java version ([#446](https://github.com/bihealth/snappy-pipeline/issues/446)) ([46b2303](https://github.com/bihealth/snappy-pipeline/commit/46b2303d4ec185409dc028fa1d1d15b3d2238507))
* fixes multiple issues with the recent cbioportal export changes ([#425](https://github.com/bihealth/snappy-pipeline/issues/425)) ([740dda5](https://github.com/bihealth/snappy-pipeline/commit/740dda551cf71e49faa31c8fd76ed46aad0dc91c))
* Fixing feature-effects compression ([#440](https://github.com/bihealth/snappy-pipeline/issues/440)) ([98825db](https://github.com/bihealth/snappy-pipeline/commit/98825db29d3583780b6543b9f5456e0be5e9978b))
* fixing gatk4_hc combine_gvcfs wrapper ([e489d93](https://github.com/bihealth/snappy-pipeline/commit/e489d93064956976af311bb6efb53d3335816fcd))
* fixing problem in running varfish-export ([#378](https://github.com/bihealth/snappy-pipeline/issues/378)) ([6d95d88](https://github.com/bihealth/snappy-pipeline/commit/6d95d885e0dda562f6812ebce65d0ee447398e33))
* gcnv contig wrapper scripts do not clean their output folders ([#443](https://github.com/bihealth/snappy-pipeline/issues/443)) ([#468](https://github.com/bihealth/snappy-pipeline/issues/468)) ([618f10d](https://github.com/bihealth/snappy-pipeline/commit/618f10d6ebdd77b40c2d65d506e685ed6139a0fe))
* gcnv wrapper regex curly-brackets parsed by python string formatting ([#489](https://github.com/bihealth/snappy-pipeline/issues/489) ) ([#501](https://github.com/bihealth/snappy-pipeline/issues/501)) ([8fb40e5](https://github.com/bihealth/snappy-pipeline/commit/8fb40e53e4f1413e9c03dc2177eb1f5008b3dbd0))
* get_output_files from ngs_mapping as reused in variant_export_external ([#346](https://github.com/bihealth/snappy-pipeline/issues/346)) ([f4fcd7e](https://github.com/bihealth/snappy-pipeline/commit/f4fcd7ea9d783c82db3d5885255645d8dc20599f))
* give Delly more memory ([#373](https://github.com/bihealth/snappy-pipeline/issues/373)) ([05f7988](https://github.com/bihealth/snappy-pipeline/commit/05f79887478e4f07b89b228893eda547a250d290))
* give melt more memory ([#372](https://github.com/bihealth/snappy-pipeline/issues/372)) ([33785b2](https://github.com/bihealth/snappy-pipeline/commit/33785b22530419fa4b628abf08918ac84c6eeb10))
* icrease allowed running time of gCNV jobs ([#349](https://github.com/bihealth/snappy-pipeline/issues/349)) ([78521f5](https://github.com/bihealth/snappy-pipeline/commit/78521f57aabdfa0876f648a01463ab4b8734ebd4))
* integration of gcnv and melt into sv_calling_wgs ([#350](https://github.com/bihealth/snappy-pipeline/issues/350)) ([5ecad1a](https://github.com/bihealth/snappy-pipeline/commit/5ecad1a230a386b818cefb6484930b122496ce0c))
* interpret skip_libraries in sv_calling ([#383](https://github.com/bihealth/snappy-pipeline/issues/383)) ([60be51d](https://github.com/bihealth/snappy-pipeline/commit/60be51da58916ff19bd80f2b4bb4277d9bbca986))
* make GNU parallel more robust with --plain and --workdir ([#328](https://github.com/bihealth/snappy-pipeline/issues/328)) ([#331](https://github.com/bihealth/snappy-pipeline/issues/331)) ([5a4c39b](https://github.com/bihealth/snappy-pipeline/commit/5a4c39b73f1e89814ba82b0293cf8c5349faa117))
* making compatible with varfish export and mehari ([#420](https://github.com/bihealth/snappy-pipeline/issues/420)) ([a40a42d](https://github.com/bihealth/snappy-pipeline/commit/a40a42da3762f8f06a89d67b46e3f817a2f135bb))
* manta germline_targeted wrapper ([#361](https://github.com/bihealth/snappy-pipeline/issues/361)) ([8e718c2](https://github.com/bihealth/snappy-pipeline/commit/8e718c23ca25c78cdc5cb344d51e879074356a1d))
* mehari wrapper vcf preprocessing does not actually fix info svlen header ([#450](https://github.com/bihealth/snappy-pipeline/issues/450)) ([#451](https://github.com/bihealth/snappy-pipeline/issues/451)) ([55146c7](https://github.com/bihealth/snappy-pipeline/commit/55146c7f50fb22f38609b40b3e91cf891bbcaf87))
* melt wrappers ([#335](https://github.com/bihealth/snappy-pipeline/issues/335)) ([66fcd46](https://github.com/bihealth/snappy-pipeline/commit/66fcd46539c0ee041788c02d7fc970ae0fe3af5f))
* missing check for `has_annotation` in somatic_variant_filtration rules ([#505](https://github.com/bihealth/snappy-pipeline/issues/505)) ([d7420f7](https://github.com/bihealth/snappy-pipeline/commit/d7420f74c5bb5863b67e591e363bf1d68d90ee7b))
* output pattern problem in gatk4_hc/genotyper wrapper ([#329](https://github.com/bihealth/snappy-pipeline/issues/329)) ([3d9555d](https://github.com/bihealth/snappy-pipeline/commit/3d9555d5cfa11bad7f2e73987ef7ed0b1fe5889c))
* paths in gcnv wrappers ([#380](https://github.com/bihealth/snappy-pipeline/issues/380)) ([fcece11](https://github.com/bihealth/snappy-pipeline/commit/fcece112fb7de1560be3f9e1b0304ee8d7014795))
* pin dependency joblib of gcnv ([#412](https://github.com/bihealth/snappy-pipeline/issues/412)) ([02fa2da](https://github.com/bihealth/snappy-pipeline/commit/02fa2daa08ee9eaf68d5091765c9fcb03b989820))
* pin dependency version in all wrapper conda env yamls ([#492](https://github.com/bihealth/snappy-pipeline/issues/492)) ([79499d9](https://github.com/bihealth/snappy-pipeline/commit/79499d9a8c16fee7ad519c8029cc94e5de31029d))
* popdel wrappers ([#362](https://github.com/bihealth/snappy-pipeline/issues/362)) ([07d73f4](https://github.com/bihealth/snappy-pipeline/commit/07d73f40bff5656ff4f27224d584e9e9142bcc41))
* re-enabling and fixing Popdel ([#359](https://github.com/bihealth/snappy-pipeline/issues/359)) ([b28ed75](https://github.com/bihealth/snappy-pipeline/commit/b28ed7531e442657480ddd446da7737c7ecdd053))
* remove "tar results" from manta wrappers ([#376](https://github.com/bihealth/snappy-pipeline/issues/376)) ([9717be7](https://github.com/bihealth/snappy-pipeline/commit/9717be72dfb0a3536378ef726cdd93750f862e9b))
* remove hard coded contig names from gcnv wrappers ([#489](https://github.com/bihealth/snappy-pipeline/issues/489) ) ([#490](https://github.com/bihealth/snappy-pipeline/issues/490)) ([4b89594](https://github.com/bihealth/snappy-pipeline/commit/4b895943f141445a22fdb6de9acf44789fdc097a))
* remove some incorrectly expected MELT output ([1e24712](https://github.com/bihealth/snappy-pipeline/commit/1e24712b2ab0cc56725ff26288998c0aa06877e2))
* rename _version.py to version.py ([#521](https://github.com/bihealth/snappy-pipeline/issues/521)) ([4fb8d8d](https://github.com/bihealth/snappy-pipeline/commit/4fb8d8d7f8371c6c4332b6ef58c5cc2a8db52225))
* set locale to C in bcftools/TMB wrapper ([#507](https://github.com/bihealth/snappy-pipeline/issues/507)) ([9789a2d](https://github.com/bihealth/snappy-pipeline/commit/9789a2d70f5f7ed025f9a66bdae03b6bb63abda7))
* setup.py had reference to removed script ([#330](https://github.com/bihealth/snappy-pipeline/issues/330)) ([30447bd](https://github.com/bihealth/snappy-pipeline/commit/30447bd0051b3c2182743ed50be0fed35719db7a))
* sleep a random time to work around GNU parallel bug ([#328](https://github.com/bihealth/snappy-pipeline/issues/328)) ([#337](https://github.com/bihealth/snappy-pipeline/issues/337)) ([de9f2df](https://github.com/bihealth/snappy-pipeline/commit/de9f2df88f2231c701f96942ba845d86e85c0d2a))
* sleep syntax in Manta wrappers ([#358](https://github.com/bihealth/snappy-pipeline/issues/358)) ([e50d4d4](https://github.com/bihealth/snappy-pipeline/commit/e50d4d49046fb5f93bbd450da4e93c61e2adeb24))
* solve issues with variant_export_external and re-enable wgs_cnv_export_external ([#354](https://github.com/bihealth/snappy-pipeline/issues/354)) ([47d9e33](https://github.com/bihealth/snappy-pipeline/commit/47d9e33033dde889f73ee244a691533fe66ba211))
* update snakemake version requirements ([#398](https://github.com/bihealth/snappy-pipeline/issues/398)) ([427803c](https://github.com/bihealth/snappy-pipeline/commit/427803c86243ca3bc4c55951020195cb548bd803))
* use local /tmp directory in gcnv wrapper ([#363](https://github.com/bihealth/snappy-pipeline/issues/363)) ([cf052c7](https://github.com/bihealth/snappy-pipeline/commit/cf052c71a81645b3a5b0d77d3be0437ea6477852))
* use only one thread for MELT ([#322](https://github.com/bihealth/snappy-pipeline/issues/322)) ([ee84aa4](https://github.com/bihealth/snappy-pipeline/commit/ee84aa4d1253f71fabbbc7ae95648aeb56c9586d))
* various fixes to make gCNV work properly for WGS and targeted ([#370](https://github.com/bihealth/snappy-pipeline/issues/370)) ([14d86e3](https://github.com/bihealth/snappy-pipeline/commit/14d86e3b312c7432e5c69d818197f54e66950833))
* write out pedigrees was missing for sv_calling_wgs ([#379](https://github.com/bihealth/snappy-pipeline/issues/379)) ([573c72f](https://github.com/bihealth/snappy-pipeline/commit/573c72fd2c99df7c533837eb9af98c6a84031594))
* zcat usage in bwa_mem2 wrapper ([#364](https://github.com/bihealth/snappy-pipeline/issues/364)) ([9cc0bb2](https://github.com/bihealth/snappy-pipeline/commit/9cc0bb2092b02f45c09ac619bd8f745cf92bf29f))


### Reverts

* not using parallel --workdir ([#328](https://github.com/bihealth/snappy-pipeline/issues/328)) ([#332](https://github.com/bihealth/snappy-pipeline/issues/332)) ([4cd1531](https://github.com/bihealth/snappy-pipeline/commit/4cd1531c17041d6c843ab55201ec54b1d2b2297c))


### Documentation

* add warning for writing new workflow steps ([#400](https://github.com/bihealth/snappy-pipeline/issues/400)) ([ac1a34f](https://github.com/bihealth/snappy-pipeline/commit/ac1a34f44addde7a654a8df991a964c2d41dda37))


### Code Refactoring

* merge wgs_{sv,cnv,mei}_calling into sv_calling_wgs ([#275](https://github.com/bihealth/snappy-pipeline/issues/275)) ([#315](https://github.com/bihealth/snappy-pipeline/issues/315)) ([fa7d0e1](https://github.com/bihealth/snappy-pipeline/commit/fa7d0e165179c31e4216082a9dd33d0a8fb06a64))
* move "bcftools roh" as report in variant_calling ([#310](https://github.com/bihealth/snappy-pipeline/issues/310)) ([#311](https://github.com/bihealth/snappy-pipeline/issues/311)) ([ee7e1fa](https://github.com/bihealth/snappy-pipeline/commit/ee7e1faffb471a800dbcc8736679f96b25f2b4f6))
* rename helper_gcnv_model_target_seq to helper_gcnv_targeted ([#303](https://github.com/bihealth/snappy-pipeline/issues/303)) ([#304](https://github.com/bihealth/snappy-pipeline/issues/304)) ([3f36328](https://github.com/bihealth/snappy-pipeline/commit/3f363282c5e62c6318aad10fd8897e76e3391d95))
* rename targeted_seq_cnv_calling to sv_calling_targeted ([#305](https://github.com/bihealth/snappy-pipeline/issues/305)) ([#306](https://github.com/bihealth/snappy-pipeline/issues/306)) ([e234160](https://github.com/bihealth/snappy-pipeline/commit/e2341602d27da33d56f864dd2c96fa3cad759e78))

## Changelog
