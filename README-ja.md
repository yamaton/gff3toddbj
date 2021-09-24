# GFF3-to-DDBJ

- [GFF3-to-DDBJ](#gff3-to-ddbj)
  * [これは何？](#これは何)
  * [セットアップ](#セットアップ)
      - [biocondaからconda環境にインストールする場合](#biocondaからconda環境にインストールする場合)
      - [GitHubソースコードからconda環境にインストールする場合](#githubソースコードからconda環境にインストールする場合)
  * [GFF3 と FASTA から DDBJ アノテーションをつくる](#gff3-と-fasta-から-ddbj-アノテーションをつくる)
    + [`gff3-to-ddbj` を動かしてみる](#gff3-to-ddbj-を動かしてみる)
  * [設定いろいろ](#設定いろいろ)
    + [メタデータファイル](#メタデータファイル)
    + [[パワーユーザ向け] Feature/Qualifier変換テーブルの編集](#パワーユーザ向け-featurequalifier変換テーブルの編集)
  * [トラブルシューティング](#トラブルシューティング)
    + [GFF3 の正当性チェック](#gff3-の正当性チェック)
    + [GFF3 と FASTA の分離（必要に応じて）](#gff3-と-fasta-の分離必要に応じて)
    + [Entry 名の正規化（必要に応じて）](#entry-名の正規化必要に応じて)
  * [当プログラムが行うこと](#当プログラムが行うこと)
  * [謝辞](#謝辞)


[TOC]

## これは何？

DDBJへの登録には指定された形式の[アノテーションファイル](https://www.ddbj.nig.ac.jp/ddbj/file-format.html#annotation)が必要です。GFF3-to-DDBJ はこの**アノテーションファイルをFASTA と GFF3ファイルから作る**プログラムです。**FASTA 単体から最小限のアノテーションファイルを作ることも可能**です。


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



## GFF3 と FASTA から DDBJ アノテーションをつくる

### `gff3-to-ddbj` を動かしてみる

まずはスクリプトを動かしてみます。

* `--gff3 <FILE>` には入力 GFF3ファイルへのパスを指定
* `--fasta <FILE>` には入力 FASTAファイルへのパスを指定
* `--metadata <FILE>` にはTOML形式のメタデータファイルを指定
* `--locus_tag_prefix <STRING>` には [BioSampleの申請時に得られたもの](https://www.ddbj.nig.ac.jp/ddbj/file-format.html#locus_tag)を指定
* `--transl_table <INT>` は [The Genetic Codes](https://www.ddbj.nig.ac.jp/ddbj/geneticcode.html) から適切な数字を選ぶ
* `--output <FILE>` には出力ファイル（＝アノテーション）のパスを指定

```shell
gff3-to-ddbj \
  --gff3 myfile.gff3 \               # この行を削除すると source, assembly_gap のみになります
  --fasta myfile.fa \                # <<必須>>
  --metadata mymetadata.toml \       # この行を削除すると作成者名などにサンプル値が入ります
  --locus_tag_prefix MYOWNPREFIX \   # この行を削除すると LOCUSTAGPREFIX_ に設定されます
  --transl_table 1 \                 # この行を削除すると 1 に設定されます
  --output myawesome_output.ann      # この行を削除すると標準出力に
```



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

* COMMON に入れる[基本情報](https://www.ddbj.nig.ac.jp/ddbj/file-format.html#annotation)

  * `SUBMITTER`, `REFERENCE` といった Features

* COMMON に入れる[メタ表記](https://www.ddbj.nig.ac.jp/ddbj/file-format.html#common)

  * 例

    ```toml
    [COMMON.assembly_gap]
    estimated_length = "unknown"
    gap_type = "within scaffold"
    linkage_evidence = "paired-ends"
    ```

  * COMMON の下に Feature を入れておくことで、DDBJ でつくられる最終記載（＝フラットファイル）に一律にQualifiersが挿入される機能があるようです。たとえば以下のようにしておくと、gff3-to-ddbjの生成するアノテーションファイルにおいて COMMON エントリの下に `assembly_gap` 以下の項目が追加されます。これは最終記載において `assembly_gap` 毎に同じ Qualifier値が挿入されることになります。


* Feature ごとに挿入する Qualifier 情報

  * 例: COMMON記法とは `[COMMON.assembly_gap]` → `[assembly_gap]` が違いだけです

    ```toml
    [assembly_gap]
    estimated_length = "unknown"
    gap_type = "within scaffold"
    linkage_evidence = "paired-ends"
    ```

    * Feature ごとのQualifier値の挿入を **gff3-to-ddbj がつくるアノテーションファイルに対して**行います。上記「COMMONに入れるメタ情報」と実質的に同じですが、アノテーションファイルの時点で Feature毎に値が入れられる点が違います。使い方は `[COMMON.assembly_gap]` ではなく `[assembly_gap]` に置き換えるだけになります。現在のところ `[source]` と `[assembly_gap]` のみ対応しています。

    * COMMON記法と併記されたばあい（たとえば `[COMMON.asembly_gap` と `[assembly_gap]` ）COMMON記法が優先されます。

さらなる例としては、アノテーションファイルとしては [WGS in COMMON](https://docs.google.com/spreadsheets/d/15gLGL5FMV8gRt46ezc2Gmb-R1NbYsIGMssB0MyHkcwE/edit#gid=1110334278) と [WGS](https://docs.google.com/spreadsheets/d/15gLGL5FMV8gRt46ezc2Gmb-R1NbYsIGMssB0MyHkcwE/edit#gid=382116224)、そして対応するメタデータファイルでは [metadata_WGS_COMMON.toml](https://github.com/yamaton/gff3toddbj/blob/main/examples/metadata/metadata_WGS_COMMON.toml) と [metadata_WGS.toml](https://github.com/yamaton/gff3toddbj/blob/main/examples/metadata/metadata_WGS.toml) などを参考にしてください。



### [パワーユーザ向け] Feature/Qualifier変換テーブルの編集

GFF3 と DDBJ アノテーションには大まかに以下のような対応があります。

1.  GFF3 の3列目 → DDBJ アノテーション 2列目 Feature
2.  GFF3 の9列目 → DDBJ アノテーション 4，5列目 Qualifier Key, Value

そしてDDBJアノテーションとして許される Features と Qualifiers の名前には規定があります [[Feature-Qualifier 一覧表](https://docs.google.com/spreadsheets/d/1qosakEKo-y9JjwUO_OFcmGCUfssxhbFAm5NXUAnT3eM/edit#gid=0)] 。

GFF3での慣習名と INSDC / DDBJ で定められた名前の橋渡しをするため、GFF3-to-DDBJは変換テーブルを用いてリネームを行っています。たとえば `five_prime_UTR` としてGFF3 の2列目に現れるものは、アノテーションファイルでは  `5'UTR`  に置き換えられます。この変換はTOMLで以下のように書きます。

```toml
[five_prime_UTR]
target = "5'UTR"
```

またQualifiersに対しては、値の前にprefixを付けることが可能です。たとえばGFF3にて `ID=foobar`となっているものを `/note` Qualifierに `ID:foobar` という値にして入れるばあいは以下のように書きます。

```toml
[ID]
target = "note"
prefix = "ID:"
```

詳しくはデフォルト設定 [translate_features.toml](https://raw.githubusercontent.com/yamaton/gff3toddbj/main/gff3toddbj/translate_features.toml) と [translate_qualifiers.toml](https://raw.githubusercontent.com/yamaton/gff3toddbj/main/gff3toddbj/translate_qualifiers.toml) をご覧ください。これらをカスタマイズしたものを使う場合には

* `--translate_features <file>` で Feature 名の変換テーブルを指定
* `--translate_qualifiers <file>` で Qualifier 名の変換テーブルを指定

のオプションがあります。呼び出しは以下のようになります。

```shell
gff3-to-ddbj \
  --gff3 myfile.gff3 \
  --fasta myfile.fa \
  --metadata mymetadata.toml \
  --locus_tag_prefix MYOWNPREFIX_ \
  --transl_table 1 \
  --translate_features translate_features.toml \       # Feature の変換テーブル指定
  --translate_qualifiers  translate_qualifiers.toml \  # Qualifier の変換テーブル指定
  --output myawesome_output.ann
```



## トラブルシューティング

### GFF3 の正当性チェック

アノテーションファイルへの変換を始めるまえに、手持ちのファイルがGFF3形式を満たしているかチェックをかけておくのが良いプラクティスです。オンラインで利用可能なものは [GFF3 online validator](http://genometools.org/cgi-bin/gff3validator.cgi) が便利です。ファイル上限が 50MB なのが玉に瑕ですが。



### GFF3 と FASTA の分離（必要に応じて）

GFF3 ファイル中に `##FASTA` を使っての塩基配列が含まれている場合には、同梱のツールを使うなどして分割してください。

```shell
split-fasta path/to/myfile.gff3 --suffix "_splitted"
```

このばあい `myfile_splitted.gff3` と `myfile_splitted.fa` の2つのファイルが作られます。



### Entry 名の正規化（必要に応じて）

DDBJ のアノテーションチェックソフトによると `=|>" []` といった文字は Entry として使えないとのことです。違反文字が含まれるときには GFF3の1列目 (= "SeqID") および FASTAのヘッダをリネームする必要があります。ステップ２でその旨のエラーを見かけたら以下を実行してください。

```shell
rename-ids \
  --gff3=path/to/foo.gff3 \     # <<必須>>
  --fasta=path/to/bar.fasta \   # <<必須>>
  --suffix="_renamed"           # この行を省略するとデフォルト値に
```

リネームの必要があるときには `foo_renamed.gff3` と `bar_renamed.fasta` の2つのファイルが作られます。無いときには `IDs are fine: No need to regularize them.` のメッセージが出るだけで終了します。



## 当プログラムが行うこと

表示形式を変えるほかに以下のようなことをしています。

* 変換テーブルに基づいた Features / Qualifiers のリネーム

* assembly_gap の検索

* /transl_table の CDS への追加

* メタデータ中の source 情報を各エントリに追加

* Qualifier 値に使われる文字の正規化

* 同じ親をもつ CDS を `join` 記法で結合

* 入力GFF3 ファイル中の mRNA と exon を `join`記法で結合し mRNA Feature とする

* CDS に 開始・終了コドンが無い場合の位置情報修正
  * 参照: [codon_start qualifier による翻訳開始の位置補正](https://www.ddbj.nig.ac.jp/ddbj/cds.html#frame)

* CDS 下の /product が値をひとつだけ持つよう変更。複数値の残りは /note へ

  * [DDBJ の/product 詳細](https://www.ddbj.nig.ac.jp/ddbj/qualifiers.html#product)を参照ください

    > * 一般名が複数ある場合でも, 複数の名称を記載しないで下さい。また, そのために不必要な区切り記号を使用しないで下さい。一般名の複数記載を希望される場合は, 代表的な名称を /product qualifier に記載し, その他の名称を /note qualifier に記載して下さい。
    > * 機能, 名称等が不明な蛋白質の場合は, hypothetical protein と記載することを推奨します。

* Qualifier値に重複があるときに冗長分を削除

* アノテーション行の並び替え

* [DDBJ の Feature-Qualifier 一覧表](https://www.ddbj.nig.ac.jp/assets/files/pdf/ddbj/fq-j.pdf)に基づいた出力情報のフィルタリング




## 謝辞
このプログラムの設計には、EMBL向けGFF3の変換ソフトである [EMBLmyGFF3](https://github.com/NBISweden/EMBLmyGFF3) のつくりを参考にさせていただきました。

