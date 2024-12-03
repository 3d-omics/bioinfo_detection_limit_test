# Changelog

## [1.5.0](https://github.com/3d-omics/mg_quant/compare/v1.4.0...v1.5.0) (2024-12-02)


### Features

* use config.yml ([586d5c9](https://github.com/3d-omics/mg_quant/commit/586d5c909f96fc40f784ae6693fec59a2cddea0c))


### Bug Fixes

* place coverm files so mag catalogue stays outside the sample name ([f540560](https://github.com/3d-omics/mg_quant/commit/f540560dbf87a42d8ffa62279bd3e8b5e5779a94))

## [1.4.0](https://github.com/3d-omics/mg_quant/compare/v1.3.0...v1.4.0) (2024-11-26)


### Features

* populate profiles ([79244a7](https://github.com/3d-omics/mg_quant/commit/79244a7aa3a49b1db99bc7c002fecd265e769924))


### Bug Fixes

* pick everything from devel ([7e98394](https://github.com/3d-omics/mg_quant/commit/7e983948aea8acb2bb512f7d6e694b127b0141d9))

## [1.3.0](https://github.com/3d-omics/mg_quant/compare/v1.2.0...v1.3.0) (2024-10-22)


### Features

* compress coverm files ([b4cd436](https://github.com/3d-omics/mg_quant/commit/b4cd4369333277def3939850a3479bd8e6d29a5d))
* cover, mags, multiqc; delete useless code; compress coverm files" ([af82953](https://github.com/3d-omics/mg_quant/commit/af829536817ed3e7feb62ae186faace0ace87a90))
* insert mg_preprocess and delete redundant code ([b1532a3](https://github.com/3d-omics/mg_quant/commit/b1532a393ad5301b1a52421ead6f3f2863a14d7b))
* pin dependencies ([711e700](https://github.com/3d-omics/mg_quant/commit/711e70058b3f25724838f3ed14ec891785b76b04))
* rename rules and add meta-rule ([ab9e99f](https://github.com/3d-omics/mg_quant/commit/ab9e99fbf050fc592557f8f9ca35d9af7a72e5f2))
* simplify coverm file names, remove singlem script ([a660727](https://github.com/3d-omics/mg_quant/commit/a660727fede2ff4fd627c90a1773a1dfe3ec0a2a))
* use helpers and preprocess modules. update params ([89915b4](https://github.com/3d-omics/mg_quant/commit/89915b416e93a33d1c55e958cf6f887c5e4d8f11))


### Bug Fixes

* add idxstats, ask for folder as output too ([de73fd3](https://github.com/3d-omics/mg_quant/commit/de73fd34b35b966df130597fe4b8f03386de99bb))
* typo ([50799a5](https://github.com/3d-omics/mg_quant/commit/50799a5762c6bfc1e93105520bbd6a92073108de))

## [1.2.0](https://github.com/3d-omics/mg_quant/compare/v1.1.2...v1.2.0) (2024-09-02)


### Features

* start adding support for nonpareil multiqc reports ([047a216](https://github.com/3d-omics/mg_quant/commit/047a21679590bbb7e122c3592c58d53e1dab591a))
* update nonpareil and multiqc so they do the jsoning and the reporting ([5fc6bdd](https://github.com/3d-omics/mg_quant/commit/5fc6bdd9f0e3d2b48dcba1aed135a774cbb6cd08))


### Bug Fixes

* add missing script ([39fbc23](https://github.com/3d-omics/mg_quant/commit/39fbc23a1dc121846dba36a139a19834b47b33f6))
* disable process substitution in fastp since multiqc can't handle it properly ([f593508](https://github.com/3d-omics/mg_quant/commit/f593508f506ae21fa04012cfd061452bd2d98e30))
* remove group from import ([6071514](https://github.com/3d-omics/mg_quant/commit/6071514aaaab119e2b6a2e749281e34f2a4128b7))
* wrong path ([e7ecb4e](https://github.com/3d-omics/mg_quant/commit/e7ecb4e71d2060fcc2370621ff6ee93338ccd26b))

## [1.1.2](https://github.com/3d-omics/mg_quant/compare/v1.1.1...v1.1.2) (2024-08-08)


### Bug Fixes

* update multiqc so it properly reports fastp ([ca1fbb2](https://github.com/3d-omics/mg_quant/commit/ca1fbb24cb0e5e95611440efde3ee971a43cff32))

## [1.1.1](https://github.com/3d-omics/mg_quant/compare/v1.1.0...v1.1.1) (2024-08-07)


### Performance Improvements

* correct kraken2 rule name in profiles ([870167c](https://github.com/3d-omics/mg_quant/commit/870167cd0ac1d93f988fdd4e4e3219cc01c6a571))

## [1.1.0](https://github.com/3d-omics/mg_quant/compare/v1.0.1...v1.1.0) (2024-08-06)


### Features

* pin dependencies ([2f53db9](https://github.com/3d-omics/mg_quant/commit/2f53db95e7e8b2bcf57eb18a1d65a2b50ea521ce))


### Bug Fixes

* add cache rule to the bottom ([2c3c9fd](https://github.com/3d-omics/mg_quant/commit/2c3c9fd76052d1365033307c127c264e38c6db96))
* nonpareil ([3ac7cc3](https://github.com/3d-omics/mg_quant/commit/3ac7cc354d63634ea9e2a4fec9665fbcc628ee3a))
* **singlem:** pin smafa version ([ef89c21](https://github.com/3d-omics/mg_quant/commit/ef89c219b1cb24f31c1da648ed90e2923632af64))

## [1.0.1](https://github.com/3d-omics/mg_quant/compare/v1.0.0...v1.0.1) (2024-06-06)


### Bug Fixes

* fix rule names in profile ([79ca5a0](https://github.com/3d-omics/mg_quant/commit/79ca5a01499f54bddf2b25d645bbbb80ff96c05e))

## 1.0.0 (2024-05-13)


### Features

* 0 threads when copying, rename group ([83b4290](https://github.com/3d-omics/mg_quant/commit/83b42901ad8676259f7c923cbc74e205b7966476))
* accelerate fastp ([e5ba830](https://github.com/3d-omics/mg_quant/commit/e5ba8306ad9c5975aed4d6e7a71899cb79b23c3f))
* **adaptors:** adaptors should be explicitely written ([cc12e74](https://github.com/3d-omics/mg_quant/commit/cc12e7486e0fbc3d36068a63451f22d517529d4b))
* add execution groups ([2d21e86](https://github.com/3d-omics/mg_quant/commit/2d21e869cd951272d6b4b8f92501e6568f59a2e2))
* add fastp to preprocess group ([df0df2e](https://github.com/3d-omics/mg_quant/commit/df0df2ec472f64597b87804c264e5ecd861aef4e))
* add groups in run_slurm ([778f2a4](https://github.com/3d-omics/mg_quant/commit/778f2a45bf4dcec4af99a95f1357ef7561b0a0f5))
* add indexes and references to cache ([16cc736](https://github.com/3d-omics/mg_quant/commit/16cc7361dcb89550ecaff075da2c89a1626e8c63))
* add micro profile ([6aad3bf](https://github.com/3d-omics/mg_quant/commit/6aad3bf76821170d6795ec3b2f46d0dcc03b0fa2))
* add profile ([7e0fc8a](https://github.com/3d-omics/mg_quant/commit/7e0fc8adab4ceec71b7df8879be2faf9a3e31af6))
* Add run_slurm [skip ci] ([a3123b1](https://github.com/3d-omics/mg_quant/commit/a3123b1293315a9cc27f6f77da9e55929101ca07))
* **bowtie2:** disable minimizer compression ([dd66c46](https://github.com/3d-omics/mg_quant/commit/dd66c464b13c948611f3f079d548437e65075e7d))
* **bowtie2:** get ram from params ([9c8162d](https://github.com/3d-omics/mg_quant/commit/9c8162de3fbb84e04165141e2138c6f74f6bf9a9))
* **bowtie2:** increase compression level ([d870e95](https://github.com/3d-omics/mg_quant/commit/d870e950c585dc6c31c12dd92c9e8fb031834b38))
* **bowtie2:** remove bottleneck ([933bece](https://github.com/3d-omics/mg_quant/commit/933becef4552f54eb05e23bdf3bd6fdc0e9ac159))
* cache bowtie2 indexes ([4ec4a55](https://github.com/3d-omics/mg_quant/commit/4ec4a55f653612ff47479c9df332d9d3aa29c517))
* cache omitting software versions ([f91122c](https://github.com/3d-omics/mg_quant/commit/f91122c36231362af7f2932cfe4e080f8fef0d6a))
* clean code ([2006df1](https://github.com/3d-omics/mg_quant/commit/2006df1ff32d0239ee57d4387404c26c294fb38c))
* consistent naming ([cb49e19](https://github.com/3d-omics/mg_quant/commit/cb49e1963c16abd6bd18b9fd48eb46a6408b3d35))
* correct column names in samples.tsv ([016f850](https://github.com/3d-omics/mg_quant/commit/016f8507be4cceaebb6f70eddb3ac3a7a44361b8))
* **coverm:** accept cram as input ([26b601b](https://github.com/3d-omics/mg_quant/commit/26b601bdc43742073c38adf862187625668789dc))
* **coverm:** refactor aggregator function ([74dcd06](https://github.com/3d-omics/mg_quant/commit/74dcd06899c21b39e97b90eab499842c6fc519df))
* cram to fastq in quantify. add defaults to profiles ([4055a4c](https://github.com/3d-omics/mg_quant/commit/4055a4c0fa899beec0d54531348b3ee46ac95e59))
* delete all the code mentioning SE reads ([57f3ac3](https://github.com/3d-omics/mg_quant/commit/57f3ac372e7701181d856b885070c5fce8577709))
* extract threads and resources to profile in preprocess ([c854f80](https://github.com/3d-omics/mg_quant/commit/c854f80982e96878003af8714ac764f370c710c2))
* extract threads and resources to profile in quantify ([5846fc6](https://github.com/3d-omics/mg_quant/commit/5846fc6ec151bc29505bb84e5909cd7accd85cea))
* extract threads and resources to profile in reference/recompress ([0292ca0](https://github.com/3d-omics/mg_quant/commit/0292ca0baadbf60ffc5581a37ffa894bc97ab88f))
* extract threads and resources to profile in report ([a46a0b5](https://github.com/3d-omics/mg_quant/commit/a46a0b5691db90bb95875b316191c74891ff7458))
* extract threads and resources to profile in samtools ([d5653ff](https://github.com/3d-omics/mg_quant/commit/d5653ffdc916ae5cb1a3048355910f4ff0c32ff5))
* **fastp:** insert extra trims to fastp ([71dee60](https://github.com/3d-omics/mg_quant/commit/71dee6086e0d1ae82ce7fbc937ac062929b04ea7))
* **fastp:** manage more elegantly the adaptors and the NAs ([0a81ad8](https://github.com/3d-omics/mg_quant/commit/0a81ad86b7bda7be1bdfe41049fde09093282a56))
* **fastp:** put an assert just in case ([c3fc4a9](https://github.com/3d-omics/mg_quant/commit/c3fc4a9e4fd8f22208d880c44d34463fc87721f7))
* **fastp:** refactor functions ([854d1e3](https://github.com/3d-omics/mg_quant/commit/854d1e3450d4fd8cceb2faadf5453de27746407d))
* **fastp:** remove temps, raise compression level ([afe88aa](https://github.com/3d-omics/mg_quant/commit/afe88aa148d1756167d8f27f2e3ce6c14641b108))
* final profiles ([ca3cde5](https://github.com/3d-omics/mg_quant/commit/ca3cde54cff3f75624db79af0e367ff6fffa5412))
* **folders:** consistent reordering ([3358c0a](https://github.com/3d-omics/mg_quant/commit/3358c0a3eb912311977d2ef66fb5cbc9f7e67013))
* **functions:** remove unused code ([4a3acfa](https://github.com/3d-omics/mg_quant/commit/4a3acfa3caa3679659f8d71f2c84da61544e9a92))
* **gha:** split into steps ([269ec70](https://github.com/3d-omics/mg_quant/commit/269ec70c639cb725128bf302ef670de5ef40fb69))
* **gha:** upload logs as artifacts on failure ([b52692a](https://github.com/3d-omics/mg_quant/commit/b52692aac518db95632c3d3db432b9a05f9e3a84))
* handle big and small bowtie2 indexes ([af2f956](https://github.com/3d-omics/mg_quant/commit/af2f9567f82d1e31bc868b7c122355a95dd8d388))
* humanize default profile ([5a4c615](https://github.com/3d-omics/mg_quant/commit/5a4c615e8cfb7a3f3847d9f12b5abd758def3c69))
* in bowtie2 map rules, clean past tmp bam files ([8a68596](https://github.com/3d-omics/mg_quant/commit/8a6859656a1a66debfbfc15fd8ed51de2929e73c))
* join heavy processes into a single group ([e52712a](https://github.com/3d-omics/mg_quant/commit/e52712afd9c0142c17103d51af371502ebbfcf13))
* kraken2 as a service ([ed02c97](https://github.com/3d-omics/mg_quant/commit/ed02c97a717de39a4fdd96e8eaaba876d50ba589))
* kraken2 as in mg_assembly ([201af9b](https://github.com/3d-omics/mg_quant/commit/201af9b2a26882b3323ac7c358ce2d4f7ecb597b))
* **kraken2:** print the hostname ([a91735f](https://github.com/3d-omics/mg_quant/commit/a91735fa768a41c6911fcfa87ce5497f5134c8a0))
* lower compression level, don't produce unpaired ([60e2c67](https://github.com/3d-omics/mg_quant/commit/60e2c6776029817b0e0301e34b65aba4d1a88ba3))
* mark singlem fastqs as temp ([62c3943](https://github.com/3d-omics/mg_quant/commit/62c3943ab29698a0cdf2221f70785f14ea7f4b8f))
* mem_mb -&gt; mem_gb ([58d085f](https://github.com/3d-omics/mg_quant/commit/58d085f8aabc001ff3f43a8f163aead64eb19e51))
* **nonpareil:** cleaner code ([345f381](https://github.com/3d-omics/mg_quant/commit/345f3818b678ced680dcf8dea9c0aff548e4f607))
* **organization:** dunder all the important files ([0c6b088](https://github.com/3d-omics/mg_quant/commit/0c6b088561d3aa9cd7ee691a43888e76de0a205c))
* **params:** remove mentions to extra fields ([f489d2f](https://github.com/3d-omics/mg_quant/commit/f489d2f8e4304f4977ad20c4797c73cd18bfe7fa))
* **pre-commit:** add R linters ([df017d8](https://github.com/3d-omics/mg_quant/commit/df017d8786c884aba4fd0e16d68ff17a8a13851a))
* **pre-commit:** add R linters ([cd26577](https://github.com/3d-omics/mg_quant/commit/cd265777bb3d6da96ba588f66b5c0acde0ae9c0e))
* **pre-commit:** raise message to the top ([aa6b764](https://github.com/3d-omics/mg_quant/commit/aa6b7649d66040a8e3ca50507c94b82bb0d52ab3))
* **preprocess:** run everything by default, and everything but singlem on testing ([5c949c2](https://github.com/3d-omics/mg_quant/commit/5c949c24083701dc29e839745f576691d34b5e56))
* raise threads for bowtie ([67ca596](https://github.com/3d-omics/mg_quant/commit/67ca596d295d2596cb3dfba30ebbda2cb260eec0))
* raise threads for bowtie ([8d02fd0](https://github.com/3d-omics/mg_quant/commit/8d02fd0b5785d3157d6661f923e58a166bf04e3a))
* **readme:** Update badges ([d5a4e1d](https://github.com/3d-omics/mg_quant/commit/d5a4e1d71220a9fa7df6714c4dea1cc67f26ae7d))
* **readme:** update with new control files ([670c818](https://github.com/3d-omics/mg_quant/commit/670c818f48ac51e1ccdf03fa42f21018a29b2e89))
* **reads:** refactor functions ([b8952ac](https://github.com/3d-omics/mg_quant/commit/b8952acec665aff1922aa600bbc9fbb87bf903fb))
* **reads:** split into subworkflows ([e05481d](https://github.com/3d-omics/mg_quant/commit/e05481df71a0f4241d3b2407211b80987a01e059))
* remove commented code. remove tmp of links ([4bc96c1](https://github.com/3d-omics/mg_quant/commit/4bc96c11c0029ad88a99c5b1f24a5298b0fbddab))
* remove dead code. cache correct files ([b1eaa94](https://github.com/3d-omics/mg_quant/commit/b1eaa94cd66697520c957a451bbae258fde56e4f))
* remove old scripts ([87bb8e0](https://github.com/3d-omics/mg_quant/commit/87bb8e01a75638dc751daa15d37b33b51170f932))
* rename and relink preprocessing ([3d73aef](https://github.com/3d-omics/mg_quant/commit/3d73aef8acc7e088a6f4c5ddeb4d1cfe9092d9a2))
* rename and relink reads ([d44c60c](https://github.com/3d-omics/mg_quant/commit/d44c60c2a8ee96aa80b886962543092a6a686957))
* rename and relink reference ([c579973](https://github.com/3d-omics/mg_quant/commit/c579973afb6f7ad1eb790963acbdfa26b163cafe))
* rename and relink report ([2e7dbbf](https://github.com/3d-omics/mg_quant/commit/2e7dbbf95c7d5609cfec48153bbb9b2107b67fc1))
* rename and relink stats ([b206ee6](https://github.com/3d-omics/mg_quant/commit/b206ee60d47da2cbc27a65762c5da5f84d4a3cd3))
* rename helpers ([374a63e](https://github.com/3d-omics/mg_quant/commit/374a63ec870d21edce991153aeeae4e63c62112f))
* rename rules and add fai to cache ([486a57f](https://github.com/3d-omics/mg_quant/commit/486a57f34052545c6da9b7494e3ead4a28910e16))
* rename rules, create bowtie2 indexes in their place, and with multiext ([9935c4b](https://github.com/3d-omics/mg_quant/commit/9935c4baa462749bdc8f0cdf3b325dc8d222a4c9))
* **reorganization:** __env__.yml =&gt; __environment__.yml ([cb685d3](https://github.com/3d-omics/mg_quant/commit/cb685d3c489db05d81b0d95edd65d3944d3c7900))
* reorganize config/params.yml ([da9af8e](https://github.com/3d-omics/mg_quant/commit/da9af8eed2588e8809f0d8fc31022d4de399d483))
* reorganize features.yml ([4e3d07c](https://github.com/3d-omics/mg_quant/commit/4e3d07c199a2bbb244fea43234046fe82eda2443))
* reorganize into preprocess and quantify ([3ce6e4c](https://github.com/3d-omics/mg_quant/commit/3ce6e4c46d741b529aa65b1694f7b3a75cae79f9))
* **reports:** double ram and use retries ([7d0e06a](https://github.com/3d-omics/mg_quant/commit/7d0e06ad1e76c5958056477810536d7041e13fe3))
* restructure helpers ([87d307d](https://github.com/3d-omics/mg_quant/commit/87d307d86023e1bf4194e616440e6a8c2dd20f3c))
* sample, library -&gt; sample_id, library_id ([59a4c5d](https://github.com/3d-omics/mg_quant/commit/59a4c5d0e4d74096762d662f2ffd4cbbb5efd4cf))
* separator now is @ ([e7d6d18](https://github.com/3d-omics/mg_quant/commit/e7d6d18cf2cd01f3721e35d729fd29390a4d132c))
* simplify kraken2 code ([0817e6a](https://github.com/3d-omics/mg_quant/commit/0817e6a485d86d0bd3fcaa7fa6cf6297e77923de))
* **singlem:** compute microbial fraction ([f274288](https://github.com/3d-omics/mg_quant/commit/f274288e7c6a15a60d9bae0020631158483b1f54))
* **singlem:** reorganize ([8793c65](https://github.com/3d-omics/mg_quant/commit/8793c65f76cc42dad3679ae5334924eeea4e5fa2))
* **singlem:** update to 0.16.0 ([5c3a698](https://github.com/3d-omics/mg_quant/commit/5c3a698b256e357166ddaea7015e73f9597bee0b))
* sort env ([6db2a90](https://github.com/3d-omics/mg_quant/commit/6db2a9043795978298236ad9c57a43f1ee9660ff))
* sort environments alphabetically ([6f2effb](https://github.com/3d-omics/mg_quant/commit/6f2effbd42d7068ffab8809c1f6cb0fb2f3d4d1c))
* sort flags in rm ([136d34e](https://github.com/3d-omics/mg_quant/commit/136d34e8d4ddc13e35964118fa14fa32d240eae8))
* test kraken2 using parallel and 1 thread ([0fd4f61](https://github.com/3d-omics/mg_quant/commit/0fd4f610ae5b8c93e1fa194e6e4b3ca43c5bb53d))
* update all environments ([18d5c28](https://github.com/3d-omics/mg_quant/commit/18d5c28095bee3a4812a75be1fdb72d705afab88))
* update envs ([0a99d8d](https://github.com/3d-omics/mg_quant/commit/0a99d8dc07955ac15c9bb4adaebfce8dee06d871))
* update minimum snakemake version ([fcb1a98](https://github.com/3d-omics/mg_quant/commit/fcb1a984bb3abf4477f99a1ed6c835f119176a65))
* use cram files instead of cram&lt;-&gt;fq conversions ([ef7dfc0](https://github.com/3d-omics/mg_quant/commit/ef7dfc0ea7e8e5308f009e381bcdca90c24ebc49))


### Bug Fixes

* **adaptors:** typo ([ffceca2](https://github.com/3d-omics/mg_quant/commit/ffceca2b67fb8771a70945d60a3fe22c2ba201d6))
* add again memory per thread in samtools ([6645c4d](https://github.com/3d-omics/mg_quant/commit/6645c4d68f7b99aecc69e0af5978c27cd4bf6c42))
* add cram to fastq rules to groups ([067240a](https://github.com/3d-omics/mg_quant/commit/067240a9496502a7e919084a61caf71aff52b3af))
* back to collating fq files ([9583a3b](https://github.com/3d-omics/mg_quant/commit/9583a3b058bc794591c45bae686788c8c07eb99b))
* back to sort -n ([a1066e6](https://github.com/3d-omics/mg_quant/commit/a1066e6a8656f6e31ea1c99c7a3016b4c34f25c1))
* back to sort by name instead of collate ([cf93e0d](https://github.com/3d-omics/mg_quant/commit/cf93e0d532d48c85450f07ff4f97f02fc7c99903))
* **bowtie2:** in sort, store the temp files next to the cram file ([bcff482](https://github.com/3d-omics/mg_quant/commit/bcff48226c7b2f1b2e65c3100c706229c4eb0b77))
* **bowtie2:** in sort, store the temp files next to the cram file ([ff2c07d](https://github.com/3d-omics/mg_quant/commit/ff2c07dbc4a639772729a958a40694688e038b08))
* cache gzi too. use samtools collate ([4b2bc6f](https://github.com/3d-omics/mg_quant/commit/4b2bc6f7c49059ed7b1fe409784de7e2448fd483))
* collating ([683a6d5](https://github.com/3d-omics/mg_quant/commit/683a6d569c500be729b2abcf477ad79bb0acdb5a))
* collating? ([7af905b](https://github.com/3d-omics/mg_quant/commit/7af905bb36a995d0ae88adff51a67a30b2117b78))
* config: right patsh , mags now with @ ([6599b97](https://github.com/3d-omics/mg_quant/commit/6599b97f698aa0e7ee2b3f517fe9b890342416c0))
* correct per library report ([d49b2fd](https://github.com/3d-omics/mg_quant/commit/d49b2fdf516d686526d9dce09cd01b295555f0de))
* correct placement for pig genome ([1c776c6](https://github.com/3d-omics/mg_quant/commit/1c776c6891c22e4a73eb86b0c1d83e1213b4cf77))
* correct rule name ([dd8056c](https://github.com/3d-omics/mg_quant/commit/dd8056ca2da805ed3fd71476a4b24609ea12c871))
* correct rule names in grouping ([9b87aee](https://github.com/3d-omics/mg_quant/commit/9b87aee32efe743721de2693e81213d2667002f4))
* **coverm:** do not filter out unaligned reads ([9e1370a](https://github.com/3d-omics/mg_quant/commit/9e1370a7cf111badf2bf6b59bb1f67c3328847a4))
* **coverm:** forgot the exclude flags ([1e9a9a1](https://github.com/3d-omics/mg_quant/commit/1e9a9a133330b5d29cfc988f4a4a4d1e888df4b2))
* do samtools sort -n so fastq are generated correctly ([0e9fd3d](https://github.com/3d-omics/mg_quant/commit/0e9fd3dea4a192b3624bac9a05c27ab885d7da69))
* extra edge in rulegraph ([4002fef](https://github.com/3d-omics/mg_quant/commit/4002fefc9fa8bbdc50f32a368674d0d158643b17))
* extract only unmapped pairs ([4030f0d](https://github.com/3d-omics/mg_quant/commit/4030f0dc6e5c90cbcec261cbd17d2ddb700fd905))
* extreme case when no magalogue is present ([777ed90](https://github.com/3d-omics/mg_quant/commit/777ed90919864ad8928d6dfa4de8f9f1a5cba0f7))
* fastp min_length, crais for coverm, more time for nonpareil ([83f0fed](https://github.com/3d-omics/mg_quant/commit/83f0fed071ba2011f0456598130344c1388dd838))
* **fastp:** raise ram to 8GB ([31f500e](https://github.com/3d-omics/mg_quant/commit/31f500ee5fd87e2d5b4ce95f8223c5be075b100b))
* **fastp:** typo ([ae2fedd](https://github.com/3d-omics/mg_quant/commit/ae2fedd4775a87fa0d3d2c6ea6a5dd63c196a7a3))
* **folders:** consistent naming ([079efe5](https://github.com/3d-omics/mg_quant/commit/079efe5f58af6c5b8eee6329ef0bc12bfbcf8f3c))
* **gha:** add again preprocess__kraken2 ([304f99d](https://github.com/3d-omics/mg_quant/commit/304f99d57e6d8f750abe95a288a3bd3ded60bd46))
* **gha:** add r packages for pre-commit ([8d6812c](https://github.com/3d-omics/mg_quant/commit/8d6812c3023098b51fd9a43d878c21e49e6523ca))
* **gha:** delete kraken2 ([d6b5a4a](https://github.com/3d-omics/mg_quant/commit/d6b5a4a0d52c6c040c988b607983d9cef1fbcbb6))
* **gha:** do not run kraken since /dev/shm is not available ([bd2937f](https://github.com/3d-omics/mg_quant/commit/bd2937fb2bb201f7c6cd44f0f50394bc4ce2e90b))
* **gha:** no comments allowed in a yaml list ([1701eca](https://github.com/3d-omics/mg_quant/commit/1701eca77d92a8d478e4ee7e998a74857132c0ac))
* **gha:** preprocessing ([29a69b7](https://github.com/3d-omics/mg_quant/commit/29a69b713d0a1431a279bd9720cbbd2d6e807175))
* **gha:** quantify, not stats ([7ffe5a0](https://github.com/3d-omics/mg_quant/commit/7ffe5a0add90daa61dd87ff2ac0dd55e15b3c6c4))
* **gha:** use custom snakemake action ([b007667](https://github.com/3d-omics/mg_quant/commit/b00766761cbac2a7e7e6f362740fc62348fae49a))
* go back to run kraken2 db copy and assignment in a single job ([dd1d550](https://github.com/3d-omics/mg_quant/commit/dd1d55082cc33049f53786bf508e0f80ba70c1a0))
* improve kraken2 ([869b223](https://github.com/3d-omics/mg_quant/commit/869b223d9e5d0961bf664b1a6e94f6190ac78b29))
* increase memory and time resources on bowtie2 and fastp [skip ci] ([8d2c09c](https://github.com/3d-omics/mg_quant/commit/8d2c09c4671a9aa8e432a87e777ba26cedd9f777))
* its mem_mb not memory ([d7c90bf](https://github.com/3d-omics/mg_quant/commit/d7c90bf2f0eba5111f1866f13051af81081f7cc5))
* kraken2 paired ([b9703ff](https://github.com/3d-omics/mg_quant/commit/b9703ff819f7b12162932d7706e655751e61ec8f))
* kraken2: back to doing it sequentially ([e855958](https://github.com/3d-omics/mg_quant/commit/e855958d54edb1e6bbbfd92ec6b087943c740479))
* **kraken2:** use whoami instead of $USER ([2ada7cf](https://github.com/3d-omics/mg_quant/commit/2ada7cfa549b287ae79a14dc25e0dd61a3aab390))
* limit fastp to 16 threads ([adbcd98](https://github.com/3d-omics/mg_quant/commit/adbcd98d003302f43596673846c5b3e03f83b713))
* linking reads should are localrules [skip ci] ([f69bce0](https://github.com/3d-omics/mg_quant/commit/f69bce03e94857131f81c2cadba34296ec16c20a))
* may give incorrect files ([9ff4ef1](https://github.com/3d-omics/mg_quant/commit/9ff4ef1c4233f1c7a1453adafa3713bf2d6ae1cc))
* **nonpareil:** f-string in params ([10b2851](https://github.com/3d-omics/mg_quant/commit/10b285148fb506278b49410379ec5f38917f6856))
* **nonpareil:** move packages between environments ([8348ee0](https://github.com/3d-omics/mg_quant/commit/8348ee08db32f2db28ea1015c41fee8692ae017a))
* **nonpareil:** wrong fq ([4a45032](https://github.com/3d-omics/mg_quant/commit/4a450320717e3e9ba54f0d7eb5b12350703d5e18))
* **nonpareil:** wrong variable names ([aedcd1e](https://github.com/3d-omics/mg_quant/commit/aedcd1e441ea04aa701fc59a2d984203c11fc27e))
* parenthesis ([8d3d352](https://github.com/3d-omics/mg_quant/commit/8d3d352c31efa104da3d2a106b6224fe53258153))
* **pre-commit:** remove R stylers. pre-commit can't find the r packages ([0119339](https://github.com/3d-omics/mg_quant/commit/0119339680229b93fddd657ea8cd9f6726f31e3c))
* **preprocess:** add csvkit to environment ([14b6c4b](https://github.com/3d-omics/mg_quant/commit/14b6c4b2a221619315fd3f2610c65c4a6e2526cb))
* proper number of threads in pigz ([c67426c](https://github.com/3d-omics/mg_quant/commit/c67426c80db217c9c20f5b53fe4d8dff39c94bc4))
* put log in the correct place ([c647e8c](https://github.com/3d-omics/mg_quant/commit/c647e8cce4c9092418ca538683ec1bfa230aeff1))
* put temp collate files next to the final cram ([ce02591](https://github.com/3d-omics/mg_quant/commit/ce02591d489d0fd73e1a5f178cf13e679e7053a4))
* raise memory in fastp fastq indexing ([1c01928](https://github.com/3d-omics/mg_quant/commit/1c019288823df19e4cf1e9864eaca4e4576d31c0))
* raise memory on bowtie2 fastq indexing and mapping [skip ci] ([eae48d9](https://github.com/3d-omics/mg_quant/commit/eae48d998fad2c61e95b988fa5d94a52c8e99fb9))
* raise memory on fastp indexing [skip ci] ([44c3fa1](https://github.com/3d-omics/mg_quant/commit/44c3fa1c302bafeaa9731687c22d1f917e490571))
* raise memory on fastq extraction for slurm ([e3e394c](https://github.com/3d-omics/mg_quant/commit/e3e394c6acbfdf49931d726912f5947e2e5c3468))
* **readme:** update and point to new path ([2f8f383](https://github.com/3d-omics/mg_quant/commit/2f8f383e7681c0c18b9c3130fc1b0ddd45ca1e77))
* **readme:** wrong config files: tsv-&gt;yml ([14de54d](https://github.com/3d-omics/mg_quant/commit/14de54d4125c598dfdf5120ee2791a2931b32048))
* remove the .fa in the renamed mags ([91cd41a](https://github.com/3d-omics/mg_quant/commit/91cd41ac98ee63c4a85d83fa2d3a368d0b9c7607))
* **report:** give 4GB to kraken2 multiqc report ([e62cbd3](https://github.com/3d-omics/mg_quant/commit/e62cbd359a70ff43a5cdf3777465590211d1e15f))
* **report:** improve readability, merge reports, remove dead code ([098ae25](https://github.com/3d-omics/mg_quant/commit/098ae25f94c22b7ea3e6a98e52b0b1b41044aa94))
* **report:** remove samtools idxstats in mag reports. They are massive ([e8cb472](https://github.com/3d-omics/mg_quant/commit/e8cb472f408a801e548fe3b8ac1d8f6dd0f58cbe))
* saner way to get unpaired from samtools [skip ci] ([585ae8f](https://github.com/3d-omics/mg_quant/commit/585ae8fff46a47316993ded617882e53ca3834ca))
* **singlem:** hardcoded output file ([434acc9](https://github.com/3d-omics/mg_quant/commit/434acc9033282e1b3932934f7e2851a3e02070dc))
* **singlem:** missing f-string ([dcdfcc5](https://github.com/3d-omics/mg_quant/commit/dcdfcc5fe0242d4bd1e25c36fdd4c66b1dbe0887))
* **singlem:** raise runtime ([a4e3951](https://github.com/3d-omics/mg_quant/commit/a4e3951857b982d07f780fd01648f3c58f915239))
* **singlem:** raise singlem pipe runtime ([eb39f00](https://github.com/3d-omics/mg_quant/commit/eb39f00dcf2153ddc363c65d0e3effa359ba56db))
* **singlem:** use 1 CPU and 16GB of RAM ([17e80d8](https://github.com/3d-omics/mg_quant/commit/17e80d8d0488b67927308f17fa2e6d4923d14b48))
* snakefmt ([7e10116](https://github.com/3d-omics/mg_quant/commit/7e10116c2d09890a7bacdb9280778c38c13e3689))
* split the collate to produce an actual bam file ([2bcd77c](https://github.com/3d-omics/mg_quant/commit/2bcd77c29f609a7bb352e84d371b264b6214e16f))
* typo ([baa7ffc](https://github.com/3d-omics/mg_quant/commit/baa7ffcbd427226da7f5a39846b906d4db1208cb))
* unfinished command ([931dfae](https://github.com/3d-omics/mg_quant/commit/931dfae112a60fea8ab8ec12e06e53bfabcbe343))
* use small index for preprocess, big index for quantify ([ce91805](https://github.com/3d-omics/mg_quant/commit/ce91805606fb809c1ea8aebfc1954ac6ff5c875f))
* wrong singlem path ([65ecfbb](https://github.com/3d-omics/mg_quant/commit/65ecfbbdc6bb1a739003430752d2af724a5d4881))


### Performance Improvements

* reduce thread usage when converting ([3d1893d](https://github.com/3d-omics/mg_quant/commit/3d1893d3d6d9a5dc82a7be63ea9bcdd89f2f5d22))
