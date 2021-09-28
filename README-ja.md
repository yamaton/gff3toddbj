# GFF3-to-DDBJ

![GitHub tag (latest by date)](https://img.shields.io/github/v/tag/yamaton/gff3toddbj?style=for-the-badge)
[![Conda (channel only)](https://img.shields.io/conda/vn/bioconda/gff3toddbj?style=for-the-badge)](https://bioconda.github.io/recipes/gff3toddbj/README.html)
[![PyPI](https://img.shields.io/pypi/v/gff3toddbj?style=for-the-badge)](https://pypi.org/project/gff3toddbj/)

* [これは何？](#これは何)
* [セットアップ](#セットアップ)
    - [biocondaからconda環境にインストールする場合](#biocondaからconda環境にインストールする場合)
    - [GitHubソースコードからconda環境にインストールする場合](#githubソースコードからconda環境にインストールする場合)
* [GFF3とFASTAからDDBJアノテーションをつくる](#gff3とfastaからddbjアノテーションをつくる)
  + [`gff3-to-ddbj` を動かしてみる](#gff3-to-ddbj-を動かしてみる)
* [当プログラムが行うこと](#当プログラムが行うこと)
* [設定いろいろ](#設定いろいろ)
  + [メタデータファイル](#メタデータファイル)
  + [[パワーユーザ向け] Features/Qualifiersのリネーム](#パワーユーザ向け-featuresqualifiersのリネーム)
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
* [謝辞](#謝辞)


[TOC]

## これは何？

DDBJへの登録には指定された形式の[アノテーションファイル](https://www.ddbj.nig.ac.jp/ddbj/file-format.html#annotation)が必要です。GFF3-to-DDBJ はこの**アノテーションファイルをFASTA と GFF3ファイルから作る**プログラムです。**FASTA 単体から最小限のアノテーションファイルを作ることも可能**です。

同種のプログラムは、NCBI登録には[GAG](https://github.com/genomeannotation/GAG)、EMBL登録には[EMBLmyGFF3](https://github.com/NBISweden/EMBLmyGFF3)などがあります。


## セットアップ

#### biocondaからconda環境にインストールする場合

```shell
# ddbjという名前でconda環境をつくってbiocondaからパッケージをインストール
conda create -n ddbj -c bioconda -c conda-forge gff3toddbj

# 環境ddbjをアクティベート
conda activate ddbj
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

# ddbjに依存パッケージ (bioconda, bcbio-gff, toml) をインストール
conda install -c bioconda -c conda-forge biopython bcbio-gff toml

# gff3-to-ddbj および付属ツールをインストール
python setup.py install
```



## GFF3とFASTAからDDBJアノテーションをつくる

### `gff3-to-ddbj` を動かしてみる

まずはスクリプトを動かしてみます。

```shell
gff3-to-ddbj \
  --gff3 myfile.gff3 \               # この行を削除すると source, assembly_gap のみになります
  --fasta myfile.fa \                # <<必須>>
  --metadata mymetadata.toml \       # この行を削除すると作成者名などにサンプル値が入ります
  --locus_tag_prefix MYOWNPREFIX \   # この行を削除すると LOCUSTAGPREFIX_ に設定されます
  --transl_table 1 \                 # この行を削除すると 1 に設定されます
  --output myawesome_output.ann      # この行を削除すると標準出力に
```

* `--gff3 <FILE>` には入力GFF3ファイルへのパスを指定
* `--fasta <FILE>` には入力FASTAファイルへのパスを指定
* `--metadata <FILE>` にはTOML形式のメタデータファイルを指定
* `--locus_tag_prefix <STRING>` には [BioSampleの申請時に得られたもの](https://www.ddbj.nig.ac.jp/ddbj/file-format.html#locus_tag)を指定
* `--transl_table <INT>` は [The Genetic Codes](https://www.ddbj.nig.ac.jp/ddbj/geneticcode.html) から適切な数字を選ぶ
* `--output <FILE>` には出力ファイル（＝アノテーション）のパスを指定



## 当プログラムが行うこと

表示形式を変えるほかに以下のようなことをしています。

* FASTAファイルの配列をSQLiteデータベースに保存（メモリ節約のため）
  * データベースファイルは処理後に削除されます

* 変換テーブルに基づいたFeatures / Qualifiersのリネーム

* assembly_gap の検索

* /transl_table のCDSへの追加

* メタデータ中のsource情報を各エントリに追加

* 同じ親をもつCDSを`join`記法で結合

* 同じ親をもつexonを`join`記法で結合し、対応するmRNAのattributesと合成してmRNA Featureとする

* CDSに開始・終了コドンが無い場合の位置情報修正
  * 参照: [codon_start qualifier による翻訳開始の位置補正](https://www.ddbj.nig.ac.jp/ddbj/cds.html#frame)

* CDS下の `/product` が値をひとつだけ持つよう変更。値が無いときには "hypothetical protein" に。複数値の残りは `/note` へ。

  * [DDBJ の/product 詳細](https://www.ddbj.nig.ac.jp/ddbj/qualifiers.html#product)を参照ください

    > * 一般名が複数ある場合でも, 複数の名称を記載しないで下さい。また, そのために不必要な区切り記号を使用しないで下さい。一般名の複数記載を希望される場合は, 代表的な名称を /product qualifier に記載し, その他の名称を /note qualifier に記載して下さい。
    > * 機能, 名称等が不明な蛋白質の場合は, hypothetical protein と記載することを推奨します。

* Qualifier値に重複があるとき冗長分を削除

* アノテーション行の並び替え

* [DDBJのFeature-Qualifier一覧表](https://www.ddbj.nig.ac.jp/assets/files/pdf/ddbj/fq-j.pdf)に基づいた出力のフィルタリング




## 設定いろいろ

### メタデータファイル

GFF3とFASTAに無い情報をアノテーションファイルに入れるためメタデータファイルを用意します。たとえば [DDBJサイトのアノテーション例](https://www.ddbj.nig.ac.jp/ddbj/file-format.html#annotation)に対応するメタデータはTOML形式で以下のように書かれます。

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

ファイルには以下の情報を入れられます。なおファイルに何も書かれていなくともプログラムは一応動作します。

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

    * FeatureごとのQualifier値の挿入を **gff3-to-ddbj がつくるアノテーションファイルに対して**行います。上記「COMMONに入れるメタ情報」と実質的に同じですが、アノテーションファイルの時点で Feature毎に値が入れられる点が違います。使い方は `[COMMON.assembly_gap]` ではなく `[assembly_gap]` に置き換えるだけになります。現在のところ `[source]` と `[assembly_gap]` のみ対応しています。

    * COMMON記法と併記されたばあい（たとえば `[COMMON.asembly_gap` と `[assembly_gap]` ）COMMON記法が優先されます。

さらなる例としては、DDBJの[サンプルアノテーション](https://www.ddbj.nig.ac.jp/ddbj/file-format.html#sample)の [WGS in COMMON](https://docs.google.com/spreadsheets/d/15gLGL5FMV8gRt46ezc2Gmb-R1NbYsIGMssB0MyHkcwE/edit#gid=1110334278) と [WGS](https://docs.google.com/spreadsheets/d/15gLGL5FMV8gRt46ezc2Gmb-R1NbYsIGMssB0MyHkcwE/edit#gid=382116224)に対応するメタデータファイル [metadata_WGS_COMMON.toml](https://github.com/yamaton/gff3toddbj/blob/main/examples/metadata/metadata_WGS_COMMON.toml) と [metadata_WGS.toml](https://github.com/yamaton/gff3toddbj/blob/main/examples/metadata/metadata_WGS.toml) を用意しています。参考にしてください。



### [パワーユーザ向け] Features/Qualifiersのリネーム

GFF3 と DDBJ アノテーションには大まかに以下のような対応があります。

1.  GFF3 の3列目 "type" → DDBJ アノテーション 2列目 Feature Key
2.  GFF3 の9列目 "attribute" → DDBJ アノテーション 4，5列目 Qualifier Key, Value

そしてDDBJアノテーションとして許されるFeaturesとQualifiersの名前には規定があります [[Feature-Qualifier 一覧表](https://docs.google.com/spreadsheets/d/1qosakEKo-y9JjwUO_OFcmGCUfssxhbFAm5NXUAnT3eM/edit#gid=0)]。

GFF3で使われる名前・値とINSDCやDDBJで定められた名前・値の橋渡しをするため、GFF3-to-DDBJはTOML形式の変換テーブルを読み込んでリネームを行っています。

以下は、デフォルトの変換テーブルで不十分なばあいのカスタマイズについてです。


#### Type / Feature keyのリネーム

GFF3 の type として現れるものがそのままFeature keyになってほしくないときには、このリネームを行います。
たとえばデフォルト設定では `five_prime_UTR` としてGFF3の3列目 ("type") に現れるものは、アノテーションファイルでは `5'UTR` というFeature keyに置き換えられます。この変換はTOMLで以下のように書きます。

```toml
[five_prime_UTR]
feature_key = "5'UTR"
```

#### Attribute / Qualifier keyのリネーム

Typeに関係なく指定のattributeについて一律リネームを行うばあいです。このケースでは値にprefixを付けることも可能です。たとえばデフォルト設定ではGFF3での `ID=foobar` というattributeはすべて `/note` qualifierとして `/note=ID:foobar` のように置き換えています。(なお Qualifier key を明示するとき `/note` のように `/`をつけて表記する慣習に従っています。ただしアノテーションではスラッシュを付けない規則のため TOMLファイル中にて `/note` のように書くことはありません。)

対応する設定はこちらです。`__ANY__` という決め打ちのtype名の下にattributeを書きます。

```toml
[__ANY__]  # この行は構造定義のため必要
[__ANY__.ID]
qualifier_key = "note"
qualifier_value_prefix = "ID:"   # qualifier_value_prefixは省略可
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

GFF3の特定の (type, attribute) を Feature key で置き換える必要に迫られることもあります。たとえばデフォルト設定では、GFF3中で `biotype="misc_RNA"` という attribute を持つ `RNA` type を `misc_RNA` Featureに変換しています。

```toml
[RNA]    # Required though redundant
[RNA.biotype]
attribute_value = "misc_RNA"
feature_key = "misc_RNA"
```

#### カスタムしたリネーム設定で `gff3-to-ddbj` を実行する

デフォルト設定 [translate_features_qualifiers.toml](https://github.com/yamaton/gff3toddbj/blob/main/gff3toddbj/translate_features_qualifiers.toml) を参考にカスタマイズしたファイルを用意してください。カスタマイズしたTOMLファイルを読み込ませるには

* `--rename_setting <FILE>`

のコマンドラインオプションを使います。これを利用した呼び出し例は以下のようになります。

```shell
gff3-to-ddbj \
  --gff3 myfile.gff3 \
  --fasta myfile.fa \
  --metadata mymetadata.toml \
  --locus_tag_prefix MYOWNPREFIX_ \
  --transl_table 1 \
  --rename_setting my_translate_features_qualifiers.toml \    # リネーム変換用ファイルを指定
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

* `--filter_setting <FILE>`

のコマンドラインオプションを使ってファイルを入力してください。


## トラブルシューティング

### GFF3の正当性チェック

アノテーションファイルへの変換を始めるまえに、手持ちのファイルがGFF3形式を満たしているかチェックをかけておくのが良いプラクティスです。オンラインで利用可能なものは [GFF3 online validator](http://genometools.org/cgi-bin/gff3validator.cgi) が便利です。ファイル上限が 50MB なのが玉に瑕ですが。



### GFF3とFASTA の分離（必要に応じて）

GFF3 ファイル中に `##FASTA` を使っての塩基配列が含まれている場合には、同梱のツールを使うなどして分割してください。

```shell
split-fasta path/to/myfile.gff3 --suffix "_splitted"
```

このばあい `myfile_splitted.gff3` と `myfile_splitted.fa` の2つのファイルが作られます。



### Entry名の正規化（必要に応じて）

DDBJ のアノテーションチェックによると `=|>" []` といった文字は Entry として使えないとのことです。違反文字が含まれるときにはアノテーションの1列目エントリ名を正規化する（＝リネームする）必要があり、同梱のツール `normalize-entry-names`が役立ちます。これはたとえば `ERS324955|SC|contig000013` というエントリ名を `ERS324955:SC:contig000013` に直します。

```shell
normalize-entry-names myannotation_output.txt
```
アノテーションファイルのエントリ名に正規化の必要があるときには `myannotation_output_renamed.txt` のファイルが作られます。無いときには `Entry names are fine: No need to normalize.` のメッセージが出て終了します。





## 謝辞
このプログラムの設計には、EMBL向けGFF3の変換ソフトである [EMBLmyGFF3](https://github.com/NBISweden/EMBLmyGFF3) のつくりを参考にさせていただきました。

