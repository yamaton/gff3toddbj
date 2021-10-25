# GFF3-to-DDBJ

![GitHub tag (latest by date)](https://img.shields.io/github/v/tag/yamaton/gff3toddbj?style=for-the-badge)
[![Conda (channel only)](https://img.shields.io/conda/vn/bioconda/gff3toddbj?style=for-the-badge)](https://bioconda.github.io/recipes/gff3toddbj/README.html)
[![PyPI](https://img.shields.io/pypi/v/gff3toddbj?style=for-the-badge)](https://pypi.org/project/gff3toddbj/)

* [これは何？](#これは何)
* [出力の「正しさ」について](#出力の正しさについて)
* [セットアップ](#セットアップ)
    - [biocondaからconda環境にインストールする場合](#biocondaからconda環境にインストールする場合)
    - [pipからconda環境にインストールする場合](#pipからconda環境にインストールする場合)
    - [GitHubソースコードからconda環境にインストールする場合](#githubソースコードからconda環境にインストールする場合)
* [GFF3とFASTAからDDBJアノテーションをつくる](#gff3とfastaからddbjアノテーションをつくる)
  + [`gff3-to-ddbj` を動かしてみる](#gff3-to-ddbj-を動かしてみる)
* [当プログラムが行うこと](#当プログラムが行うこと)
* [設定いろいろ](#設定いろいろ)
  + [メタデータファイル](#メタデータファイル)
  + [[パワーユーザ向け] FeaturesとQualifiersのリネーム](#パワーユーザ向け-featuresとqualifiersのリネーム)
    + [Type / Feature keyのリネーム](#type--feature-keyのリネーム)
    + [Attribute / Qualifier keyのリネーム](#attribute--qualifier-keyのリネーム)
    + [指定のTypeをQualifier付きFeatureに置き換え](#指定のtypeをqualifier付きfeatureに置き換え)
    + [指定の (type, attribute) をFeatureに置き換え](#指定の-type-attribute-をfeatureに置き換え)
    + [カスタムしたリネーム設定で `gff3-to-ddbj` を実行する](#カスタムしたリネーム設定で-gff3-to-ddbj-を実行する)
  + [[パワーユーザ向け] Features-Qualifiers出力のカスタマイズ](#パワーユーザ向け-features-qualifiers出力のカスタマイズ)
* [トラブルシューティング](#トラブルシューティング)
  + [GFF3の正当性チェック](#gff3の正当性チェック)
  + [GFF3とFASTAの分離（必要に応じて）](#gff3とfastaの分離必要に応じて)
  + [Entry名の正規化（必要に応じて）](#entry名の正規化必要に応じて)
* [既知の問題](#既知の問題)
* [謝辞](#謝辞)


[TOC]

## これは何？

DDBJへの登録には指定された形式の[アノテーションファイル](https://www.ddbj.nig.ac.jp/ddbj/file-format.html#annotation)が必要です。GFF3-to-DDBJ はこの**アノテーションファイルをFASTA と GFF3ファイルから作る**プログラムです。**FASTA 単体から最小限のアノテーションファイルを作ることも可能**です。

同種のプログラムは、NCBI登録には[table2asn](https://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/table2asn_GFF/)や[GAG](https://github.com/genomeannotation/GAG)、EMBL登録には[EMBLmyGFF3](https://github.com/NBISweden/EMBLmyGFF3)などがあります。

[テストフォルダ](https://github.com/yamaton/gff3toddbj/tree/main/tests/golden)に入出力のファイルの例があるのでご覧ください。.annの拡張子付きがDDBJアノテーションファイルになります。同じファイル名のFASTAとGFF3から生成しています。


## 出力の「正しさ」について

DDBJアノテーション出力には幾多の守るべきルールはあるものの、正解には曖昧さが伴います。そのため当ツールではRefSeqの出すGFF3とGenBank形式の対応をもって「正解」と定義しています。（ただし指針がDDBJと異なるばあいにはDDBJのほうを優先しています。）これに基づいて、RefSeqのGenBank形式をできるだけシンプルにDDBJアノテーション形式に変換した正解例を、当ツールで作られるDDBJアノテーション出力とを比べることで評価をしています。[評価の詳細ページはこちら](https://github.com/yamaton/gff3toddbj/tree/main/evaluation)。

なお評価の副産物としてのGenBankからDDBJアノテーションの変換ツール `genbank-to-ddbj` を同梱しています。良い品質のGenBank形式があるときにはGFF3から作るよりも有用かもしれません。

評価の一環として、DDBJ公開の[Parser](https://www.ddbj.nig.ac.jp/ddbj/parser.html)も利用しています。


## セットアップ

#### biocondaからconda環境にインストールする場合

```shell
# ddbjという名前でconda環境をつくってbiocondaからパッケージをインストール
conda create -n ddbj -c bioconda -c conda-forge gff3toddbj

# 環境ddbjをアクティベート
conda activate ddbj
```

### pipからconda環境にインストールする場合

```shell
# ddbjという名前でconda環境をつくってpipコマンドをインストール
conda create -n ddbj pip

# 環境ddbjをアクティベート
conda activate ddbj

# samtoolsに含まれる実行ファイルbgzipが必要
conda install -c bioconda samtools

# pipからgff3toddbjをインストール
pip install gff3toddbj
```

#### GitHubソースコードからconda環境にインストールする場合

```shell
# ソースをzipでダウンロード
wget https://github.com/yamaton/gff3_to_ddbj/archive/refs/heads/main.zip

# zipを展開、リネーム
unzip main.zip && mv gff3toddbj-main gff3toddbj && cd gff3toddbj

# ddbjという名前でconda環境をつくる
conda create -n ddbj

# 環境ddbjをアクティベート
conda activate ddbj

# ddbjに依存パッケージ (bioconda, bcbio-gff, toml, setuptools) をインストール
conda install -c bioconda -c conda-forge biopython bcbio-gff toml pysam samtools pip build

# gff3-to-ddbj および付属ツールをインストール
python -m build && pip install -e ./
```



## GFF3とFASTAからDDBJアノテーションをつくる

### `gff3-to-ddbj` を動かしてみる

まずはスクリプトを動かしてみます。

```shell
gff3-to-ddbj \
  --gff3 myfile.gff3 \               # この行を削除すると source, assembly_gap のみになります
  --fasta myfile.fa \                # <<必須>>
  --metadata mymetadata.toml \       # この行を削除するとCOMMON項目無しのデフォルト値が入ります
  --locus_tag_prefix MYOWNPREFIX \   # この行を削除すると LOCUSTAGPREFIX_ に設定されます
  --transl_table 1 \                 # この行を削除すると 1 に設定されます
  --output myawesome_output.ann      # この行を削除すると標準出力に
```

* `--gff3 <FILE>` には入力GFF3ファイルへのパスを指定
* `--fasta <FILE>` には入力FASTAファイルへのパスを指定
* `--metadata <FILE>` にはTOML形式の[メタデータファイル](#メタデータファイル)を指定
* `--locus_tag_prefix <STRING>` には [BioSampleの申請時に得られたもの](https://www.ddbj.nig.ac.jp/ddbj/file-format.html#locus_tag)を指定
* `--transl_table <INT>` は [The Genetic Codes](https://www.ddbj.nig.ac.jp/ddbj/geneticcode.html) から適切な数字を選ぶ
* `--output <FILE>` には出力ファイル（＝アノテーション）のパスを指定



## 当プログラムが行うこと

表示形式を変えるほかに以下のようなことをしています。

* FASTAファイルがgzipで圧縮されているばあい[bgzip](https://www.htslib.org/doc/bgzip.html)圧縮を作成
  * FASTAのインデックス化とメモリ節約のため
  * `myname_bgzip.fa.gz` のように `_bgzip` の付いたファイルが作成されます
  * bgzipfファイルはgzipと互換性があるのでそのまま他での利用も可能です

* Features / Qualifiersの[変換テーブル](https://github.com/yamaton/gff3toddbj/blob/main/gff3toddbj/translate_features_qualifiers.toml)に基づいたリネーム
  * **これが当プログラムのコア機能になります**
  * [Sequence Ontology （略称SO）](http://sequenceontology.org)をベースに、実際例により拡張したものとなっています
  * たとえば GFF3 の3列目("type")に現れる "transcript" は ([SO:0000673](http://sequenceontology.org/browser/current_svn/term/SO:0000673))にて "INSDC_feature:misc_RNA" とあるため `misc_RNA` feature に置き換えられます。

* `assembly_gap` の検索・追加

* `/transl_table` の`CDS`への追加

* [メタデータファイル](#メタデータファイル)中の`source`情報を各エントリに追加

* GFF3が `Is_circular=true` を含むばあい `TOPOLOGY` を `/circular` 付で追加して環状ゲノムであることを明示
  * このとき[始点・終点をまたぐ featureのlocation処理](https://www.ddbj.nig.ac.jp/faq/ja/how-to-describe-location-circular-genome.html)も行う

* 同じ親をもつfeatureを`join`記法で結合
  *  `CDS`, `exon`, `mat_peptide`, `V_segment`, `C_region`, `D-loop`, `misc_feature` に対して適用
  * 例外として `gene`を直接の親に持つ `exon` は結合しない

* RNAの位置情報を、ぶら下がる結合されたexonの位置情報で置き換える

* CDSに開始・終了コドンが無い場合マークアップ (`<`, `>`) で位置の補正
  * 参照: [codon_start qualifier による翻訳開始の位置補正](https://www.ddbj.nig.ac.jp/ddbj/cds.html#frame)

* CDS下の `/product` が値をひとつだけ持つよう変更。値が無いときには "hypothetical protein" に。複数値の残りは `/note` へ。

  * 参照: [DDBJ の/product 詳細](https://www.ddbj.nig.ac.jp/ddbj/qualifiers.html#product)

    > * 一般名が複数ある場合でも, 複数の名称を記載しないで下さい。また, そのために不必要な区切り記号を使用しないで下さい。一般名の複数記載を希望される場合は, 代表的な名称を /product qualifier に記載し, その他の名称を /note qualifier に記載して下さい。
    > * 機能, 名称等が不明な蛋白質の場合は, hypothetical protein と記載することを推奨します。

* `gene` featureが `/gene` または `/gene_synonym` を持つ場合、同qualifiersを下層featuresに対しコピー。

* `/gene` が値を複数持つときには、ひとつだけを `/gene` の値とする。残りを `/gene_synonym` の値に。
  * 参照: [Qualifier Key の定義: /gene](https://www.ddbj.nig.ac.jp/ddbj/qualifiers.html#gene).

* Qualifier値に重複があるとき冗長分を削除

* アノテーション行の並び替え
  * (開始位置、Feature Keyによる優先度、終了位置) をもってソートします
  * 優先度はこちらで[定義](https://github.com/yamaton/gff3toddbj/blob/1cea725cca2a8f3edb45bac45d7983e255285d5e/gff3toddbj/transforms.py#L763)。`source`, `TOPOLOGY` が先頭にくるようにしています。

* [DDBJのFeature-Qualifier一覧表](https://www.ddbj.nig.ac.jp/assets/files/pdf/ddbj/fq-j.pdf)に基づいた出力のフィルタリング
  * 多くのGFF3に記載されている `gene` feature はここで無くなります。注意。
  * フィルタによって捨てられる項目は、実行時に標準エラー出力として以下のようにずらっと並びます。
    ```
    WARNING: [Discarded] feature ------->  gene  <-------    (count: 49911)
    WARNING: [Discarded] feature ------->  cDNA_match  <-------      (count: 10692)
    WARNING: [Discarded] feature ------->  match  <-------   (count: 101)
    WARNING: [Discarded] feature ------->  sequence_conflict  <-------   (count: 81)
    WARNING: [Discarded] (Feature, Qualifier) = (source, db_xref)    (count: 687)
    WARNING: [Discarded] (Feature, Qualifier) = (source, Name)   (count: 687)
    WARNING: [Discarded] (Feature, Qualifier) = (source, gbkey)      (count: 687)
    WARNING: [Discarded] (Feature, Qualifier) = (source, genome)     (count: 685)
    WARNING: [Discarded] (Feature, Qualifier) = (mRNA, Parent)   (count: 57304)
    WARNING: [Discarded] (Feature, Qualifier) = (mRNA, db_xref)      (count: 114608)
    ```




## 設定いろいろ

### メタデータファイル

GFF3とFASTAに無い情報をアノテーションファイルに入れるためメタデータファイルを用意します。たとえば [DDBJサイトのアノテーション例](https://www.ddbj.nig.ac.jp/ddbj/file-format.html#annotation)に対応するメタデータはTOML形式で以下のように書かれます。なお省略なしのサンプルは[こちら](https://github.com/yamaton/gff3toddbj/blob/main/examples/metadata/metadata_ddbj_example.toml)になります。

```toml
[COMMON]
[COMMON.SUBMITTER]
ab_name = [
    "Robertson,G.R.",
    "Mishima,H."
]
contact = "Hanako Mishima"
email = "mishima@ddbj.nig.ac.jp"

## ... 中略 ...

[COMMON.COMMENT]
line = [
    "Please visit our website URL",
    "http://www.ddbj.nig.ac.jp/"
]

```

メタデータファイルには以下の情報を入れられます。なおファイルが空でもプログラムは一応動作します。

* COMMONに入れる[基本情報](https://www.ddbj.nig.ac.jp/ddbj/file-format.html#annotation)

  * `SUBMITTER`, `REFERENCE` といった Features

* COMMONに入れる[メタ表記](https://www.ddbj.nig.ac.jp/ddbj/file-format.html#common)

  * 例

    ```toml
    [COMMON.assembly_gap]
    estimated_length = "unknown"
    gap_type = "within scaffold"
    linkage_evidence = "paired-ends"
    ```

  * COMMONの下にFeatureを入れておくことで、DDBJでつくられる最終記載（＝フラットファイル）に一律にQualifiersが挿入される機能があるようです。たとえば以下のようにしておくと、gff3-to-ddbjの生成するアノテーションファイルにおいてCOMMONエントリの下に `assembly_gap` 以下の項目が追加されます。これは最終記載において `assembly_gap` 毎に同じQualifier値が挿入されることになります。


* Featureごとに挿入するQualifier情報

  * 例: COMMON記法とは `[COMMON.assembly_gap]` → `[assembly_gap]` が違うだけです

    ```toml
    [assembly_gap]
    estimated_length = "unknown"
    gap_type = "within scaffold"
    linkage_evidence = "paired-ends"
    ```

    * FeatureごとのQualifier値の挿入を **gff3-to-ddbj がつくるアノテーションファイルに対して**行います。現在のところ `[source]` と `[assembly_gap]` にのみ対応しています。上記「COMMONに入れるメタ情報」と実質的に同じことですが、アノテーションファイルの時点でFeature毎に値が入る点が違います。使い方は `[COMMON.assembly_gap]` を `[assembly_gap]` に置き換えるだけになります。

`gff3-to-ddbj` 実行時にメタデータファイルが指定されない場合は、「とりあえず」として用意された[デフォルト設定](https://github.com/yamaton/gff3toddbj/blob/main/gff3toddbj/metadata_without_COMMON.toml)がロードされます。

さらなる例としては、DDBJの[サンプルアノテーション](https://www.ddbj.nig.ac.jp/ddbj/file-format.html#sample)の [WGS in COMMON](https://docs.google.com/spreadsheets/d/15gLGL5FMV8gRt46ezc2Gmb-R1NbYsIGMssB0MyHkcwE/edit#gid=1110334278) と [WGS](https://docs.google.com/spreadsheets/d/15gLGL5FMV8gRt46ezc2Gmb-R1NbYsIGMssB0MyHkcwE/edit#gid=382116224)に対応するメタデータファイル [metadata_WGS_COMMON.toml](https://github.com/yamaton/gff3toddbj/blob/main/examples/metadata/metadata_WGS_COMMON.toml) と [metadata_WGS.toml](https://github.com/yamaton/gff3toddbj/blob/main/examples/metadata/metadata_WGS.toml) を用意しています。参考にしてください。



### [パワーユーザ向け] FeaturesとQualifiersのリネーム

GFF3 と DDBJ アノテーションには大まかに以下のような対応があります。

1.  GFF3 の3列目 "type" → DDBJ アノテーション 2列目 Feature Key
2.  GFF3 の9列目 "attribute" → DDBJ アノテーション 4，5列目 Qualifier Key, Value

そしてDDBJアノテーションとして許されるFeaturesとQualifiersの名前には規定があります [[Feature-Qualifier 一覧表](https://docs.google.com/spreadsheets/d/1qosakEKo-y9JjwUO_OFcmGCUfssxhbFAm5NXUAnT3eM/edit#gid=0)]。

GFF3で使われる名前・値とINSDCやDDBJで定められた名前・値の橋渡しをするため、GFF3-to-DDBJはTOML形式の変換設定 を読み込んでリネームを行っています。[デフォルトの変換設定](https://github.com/yamaton/gff3toddbj/blob/main/gff3toddbj/translate_features_qualifiers.toml)は[the Sequence Ontology (SO)](http://sequenceontology.org/browser)や類似ツールを参考に頑張って作ってありますが、まだまだ改善の余地があるでしょうし、またGFF3アノテーションのSequence Ontology準拠レベルにも依存すると思われます。

以下は、デフォルトの変換テーブルで不十分な場合のカスタマイズについてです。


#### Type / Feature keyのリネーム

GFF3 の type として現れるものがそのままFeature keyになってほしくないときには、このリネームを行います。
たとえばデフォルト設定では `five_prime_UTR` としてGFF3の3列目 ("type") に現れるものは、アノテーションファイルでは `5'UTR` というFeature keyに置き換えられます。この変換はTOMLで以下のように書きます。

```toml
[five_prime_UTR]
feature_key = "5'UTR"
```

#### Attribute / Qualifier のリネーム

Typeに関係なく指定のattributeについて一律リネームを行うばあいです。このケースでは値にprefixを付けることも可能です。たとえばデフォルト設定ではGFF3での `ID=foobar` というattributeはすべて `/note` qualifierとして `/note=ID:foobar` のように置き換えています。(なお Qualifier key を明示するとき `/note` のように `/`をつけて表記する慣習に従っています。ただしアノテーションではスラッシュを付けない規則のため TOMLファイル中にて `/note` のように書くことはありません。)

対応する設定はこちらです。`__ANY__` という決め打ちのtype名の下にattributeを書きます。

```toml
[__ANY__.ID]
qualifier_key = "note"
qualifier_value_prefix = "ID:"   # qualifier_value_prefixは省略可
```

またqualifierキーと値の両方を設定するには以下のように書きます。たとえば `/pseudo` DDBJ新規登録において非推奨のため、 `/pseudo=true`といった値に関係なく `/pseudogene="unknown"` で置き換える措置をデフォルトで行っています。

```toml
# /pseudo は常に /pseudogene="unknown" に置き換える
[__ANY__.pseudo]
qualifier_key = "pseudogene"
qualifier_value = "unknown"
```


#### 指定のTypeをQualifier付きFeatureに置き換え

GFF3の特定typeについて、Qualifier key, valueの付いたFeatureで置き換えることもできます。たとえばデフォルト設定では `snRNA` というtypeを `/ncRNA_class="snRNA"` というqualifierを付けたFeature `ncRNA` に変換しています。

```toml
[snRNA]
feature_key = "ncRNA"
qualifier_key = "ncRNA_class"
qualifier_value = "snRNA"
```


#### 指定の (type, attribute) をFeatureに置き換え

GFF3にてとある (type名, attributeキー, attribute value) があるとき、これを Feature key として表す必要に迫られることがあります。たとえばデフォルト設定では、GFF3の３列目が `RNA` かつ９列目が `biotype="misc_RNA"` であるとき `misc_RNA` Featureとして扱うよう変換をかけています。マッチさせたい単語をドット区切りで、featureキー、qualifierキー、qualifier値の順に書きます。

```toml
[RNA.biotype.misc_RNA]
feature_key = "misc_RNA"
```

#### カスタムしたリネーム設定で `gff3-to-ddbj` を実行する

デフォルト設定 [translate_features_qualifiers.toml](https://github.com/yamaton/gff3toddbj/blob/main/gff3toddbj/translate_features_qualifiers.toml) を参考にカスタマイズしたファイルを用意してください。カスタマイズしたTOMLファイルを読み込ませるには

* `--config_rename <FILE>`

のコマンドラインオプションを使います。これを利用した呼び出し例は以下のようになります。

```shell
gff3-to-ddbj \
  --gff3 myfile.gff3 \
  --fasta myfile.fa \
  --metadata mymetadata.toml \
  --locus_tag_prefix MYOWNPREFIX_ \
  --transl_table 1 \
  --config_rename my_translate_features_qualifiers.toml \    # リネーム変換用ファイルを指定
  --output myawesome_output.ann
```


### [パワーユーザ向け] Features-Qualifiers出力のカスタマイズ

DDBJの登録にて推奨される(Feature, Qualifier)の組は限定されているため、アノテーションの出力にフィルタリングを行っています。
デフォルトでは[DDBJのFeature/Qualifier 対応一覧表](https://www.ddbj.nig.ac.jp/assets/files/pdf/ddbj/fq-j.pdf)に準じた[TOML形式の設定](https://github.com/yamaton/gff3toddbj/blob/main/gff3toddbj/ddbj_filter.toml)を読みこみます。このファイルでは以下のような構造を取ります。

```toml
CDS = [
"EC_number",
"inference",
"locus_tag",
"note",
"product",
]

exon = [
"gene",
"locus_tag",
"note",
]
```

出力として許される Features および Qualifiers を、左辺にFeature名、右辺にQualifier名を要素にもつリストを書くことで表しています。

たとえば上記の設定では、出力は `CDS` および `exon` のFeatures のみ、そして対応する Qualifiers は記述されたキーに限定されます。
この機能をカスタマイズするばあいには、TOMLファイルを編集したうえで

* `--config_filter <FILE>`

のコマンドラインオプションを使ってファイルを入力してください。


## トラブルシューティング

### GFF3の正当性チェック

アノテーションファイルへの変換を始めるまえに、手持ちのファイルがGFF3形式を満たしているかチェックをかけておくのが良いプラクティスです。オンラインで利用可能なものは [GFF3 online validator](http://genometools.org/cgi-bin/gff3validator.cgi) が便利です。ファイル上限が 50MB なのが玉に瑕です。



### GFF3とFASTA の分離

GFF3 ファイル中に `##FASTA` を使っての塩基配列が含まれている場合にはその旨のエラーが出ます。同梱のツールを使うなどして分割してください。

```shell
split-fasta path/to/myfile.gff3 --suffix "_splitted"
```

このばあい `myfile_splitted.gff3` と `myfile_splitted.fa` の2つのファイルが作られます。



### Entry名の正規化

DDBJ のアノテーションチェックによると `=|>" []` といった文字は Entry として使えないとのことです。違反文字が含まれるときにはアノテーションの1列目エントリ名を正規化する（＝リネームする）必要があり、同梱のツール `normalize-entry-names`が役立ちます。これはたとえば `ERS324955|SC|contig000013` というエントリ名を `ERS324955:SC:contig000013` に直します。

```shell
normalize-entry-names myannotation_output.ann
```
アノテーションファイルのエントリ名に正規化の必要があるときには `myannotation_output_renamed.ann` のファイルが作られます。無いときには `Entry names are fine: No need to normalize.` のメッセージが出て終了します。


## 既知の問題

* `/trans_splicing` があるときの位置補正およびFeature `join()`
* `/transl_except` が端に来るときの位置補正
* `/exception` があるときの `/translation` 対応
* 位置表記 `123^124` がフラットファイルに記載されるようなケースのGFF3側での表現・処理
* 現段階は仕様を固めているところなので、処理速度が後回しになっています


## 謝辞
このプログラムの設計には、EMBL向けGFF3の変換ソフトである [EMBLmyGFF3](https://github.com/NBISweden/EMBLmyGFF3) のつくりを参考にさせていただきました。

