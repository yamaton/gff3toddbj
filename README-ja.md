# GFF3-to-DDBJ

[TOC]

## これは何？

DDBJへの登録には指定された形式の[アノテーションファイル](https://www.ddbj.nig.ac.jp/ddbj/file-format.html#annotation)が必要です。GFF3-to-DDBJ は **このアノテーションファイルをFASTA と GFF3ファイルから作る**プログラムです。**FASTA 単体から最小限のアノテーションを作ることも可能**です。


## セットアップ

#### [審査待ち] biocondaからconda環境にインストールする場合

```shell
# ddbjという名前でconda環境をつくってbiocondaからパッケージをインストール
# 現在 (2021-09-14) 審査待ちのためごちゃついてます
$ conda create -n ddbj -c bioconda -c conda-forge -c https://168588-42372094-gh.circle-artifacts.com/0/tmp/artifacts/packages gff3toddbj

# 環境ddbjをアクティベート
$ conda activate ddbj
```



#### GitHubソースコードからconda環境にインストールする場合

```shell
# ソースをzipでダウンロード
$ wget https://github.com/yamaton/gff3_to_ddbj/archive/refs/heads/main.zip

# zipを展開、リネーム
$ unzip main.zip && mv gff3toddbj-main gff3toddbj && cd gff3toddbj

# ddbjという名前でconda環境をつくる
$ conda create -n ddbj

# 環境ddbjをアクティベート
$ conda activate ddbj

# ddbjに依存パッケージ (bioconda, bcbio-gff, toml) をインストール
$ conda install -c bioconda -c conda-forge biopython bcbio-gff toml

# gff3-to-ddbj および付属ツールをインストール
$ python setup.py install
```



## GFF3 と FASTA から DDBJ アノテーションをつくる

### `gff3-to-ddbj` を動かしてみる

まずはスクリプトを動かしてみます。

* `--gff3 <FILE>` には入力 GFF3ファイルへのパスを指定
* `--fasta <FILE>` には入力 FASTAファイルへのパスを指定
* `--config <FILE>` にはTOMLでの設定ファイルを指定
* `--locus_tag_prefix <STRING>` には [BioSampleの申請時に得られたもの](https://www.ddbj.nig.ac.jp/ddbj/file-format.html#locus_tag)を指定
* `--transl_table <INT>` は [The Genetic Codes](https://www.ddbj.nig.ac.jp/ddbj/geneticcode.html) から適切な数字を選ぶ
* `--output <FILE>` には出力ファイル（＝アノテーション）のパスを指定

```shell
gff3-to-ddbj \
  --gff3 myfile.gff3 \               # この行を削除するとそれなりの出力に
  --fasta myfile.fa \                # <<必須>>
  --config myconfig.toml \           # この行を削除するとそれなりの出力に
  --locus_tag_prefix MYOWNPREFIX \   # この行を削除すると LOCUSTAGPREFIX_ に設定
  --transl_table 1 \                 # この行を削除すると 1 と設定
  --output myawesome_output.ann      # この行を削除すると標準出力に
```



## 設定いろいろ

### 設定ファイル

GFF3とFASTAに無い情報を追加するためTOML設定ファイルの用意をお勧めします。たとえば [DDBJのアノテーション例](https://www.ddbj.nig.ac.jp/ddbj/file-format.html#annotation)に対応する設定は[このように](https://github.com/yamaton/gff3toddbj/blob/main/examples/configs/config_ddbj_example.toml)書けます。

設定ファイルには以下の情報をTOML形式で入れていただきます。が、**省略してもプログラムは動作します。**

* COMMON に入れる[基本情報](https://www.ddbj.nig.ac.jp/ddbj/file-format.html#annotation)

* COMMON に入れる[メタ表記](https://www.ddbj.nig.ac.jp/ddbj/file-format.html#common)

  * COMMON の下に Feature を入れておくことで、DDBJ でつくられるフラットファイルに自動的に一律に値が挿入される機能があるそうです。たとえば以下のようにしておくと `assembly_gap` ごとに同じ Qualifier値が挿入されることになります。

    ```toml
    [COMMON.assembly_gap]
    estimated_length = "unknown"
    gap_type = "within scaffold"
    linkage_evidence = "paired-ends"
    ```

* Feature ごとに挿入する Qualifier 情報

  * Feature ごとのQualifier値の挿入を **gff3-to-ddbj がつくるアノテーションファイルに対して**行います。上記「COMMONに入れるメタ情報」と実質的に同じですが、メタ表記に対するDDBJ側での変換作業が未確認のためこの機能を入れています。使い方は`[COMMON.assembly_gap]` ではなく  `[assembly_gap]` に置き換えるだけです。

    ```toml
    [assembly_gap]
    estimated_length = "unknown"
    gap_type = "within scaffold"
    linkage_evidence = "paired-ends"
    ```

アノテーションにおける、COMMON下のメタ表記とFeatureごとの挿入の比較は [EST in COMMON](https://docs.google.com/spreadsheets/d/15gLGL5FMV8gRt46ezc2Gmb-R1NbYsIGMssB0MyHkcwE/edit#gid=633379952) vs [EST](https://docs.google.com/spreadsheets/d/15gLGL5FMV8gRt46ezc2Gmb-R1NbYsIGMssB0MyHkcwE/edit#gid=1753678626)の例、そして対応する設定ファイル (config_EST_in_COMMON.toml vs config_EST.toml)  が分かりやすいです。



### [パワーユーザ向け] Feature/Qualifier変換テーブルの編集

GFF3 と DDBJ アノテーションには大まかに以下のような対応があります。

1.  GFF3 の3列目 $\to$ DDBJ アノテーション 2列目 Feature
2.  GFF3 の9列目 $\to$ DDBJ アノテーション 4，5列目 Qualifier Key, Value

いっぽうでDDBJアノテーションとして許される Features と Qualifiers の名前にには規定があります。（参照: Feature-Qualifier 一覧表。）GFF3における慣習と INSDC / DDBJ で定められた名前の橋渡しをするため、変換テーブルを元に `gff3-to-ddbj`は Features と Qualifiers においてそれぞれリネームを行っています。たとえば `five_prime_UTR` という名前で GFF3 の2列目に現れるものは、 `5'UTR`という Feature に置き換えられます。

デフォルトの変換テーブルは [translate_features.toml](https://raw.githubusercontent.com/yamaton/gff3toddbj/main/gff3toddbj/translate_features.toml) と [translate_qualifiers.toml](https://raw.githubusercontent.com/yamaton/gff3toddbj/main/gff3toddbj/translate_qualifiers.toml) にて定めています。これらをカスタマイズして使うばあいには

* `--translate_features <file>` で Features の変換テーブルを指定
* `--translate_qualifiers <file>` で Qualifiers の変換テーブルを指定

してください。呼び出しは以下のようになります。

```shell
gff3-to-ddbj \
  --gff3 myfile.gff3 \
  --fasta myfile.fa \
  --config myconfig.toml \
  --locus_tag_prefix MYOWNPREFIX_ \
  --transl_table 1 \
  --translate_features translate_features.toml \      # Feature Keys の変換テーブルを指定
  --translate_qualifiers  translate_qualifiers.toml \ # Qualifier Keys の変換テーブルを指定
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

* Start codon まわり整合性チェック（ いまのところ /codon_start=1 で**ない**時のみ）

* CDS 下の /product が値をひとつだけ持つよう変更。複数値の残りは /note へ

  * [DDBJ の/product 詳細](https://www.ddbj.nig.ac.jp/ddbj/qualifiers.html#product)を参照ください

    > * 一般名が複数ある場合でも, 複数の名称を記載しないで下さい。また, そのために不必要な区切り記号を使用しないで下さい。一般名の複数記載を希望される場合は, 代表的な名称を /product qualifier に記載し, その他の名称を /note qualifier に記載して下さい。
    > * 機能, 名称等が不明な蛋白質の場合は, hypothetical protein と記載することを推奨します。

* Qualifier値に重複があるときに冗長分を削除

* アノテーション行の並び替え

* [DDBJ の Feature-Qualifier 一覧表](https://www.ddbj.nig.ac.jp/assets/files/pdf/ddbj/fq-j.pdf)に基づいた出力情報のフィルタリング




## 謝辞
このプログラムの設計には、EMBL向けGFF3の変換ソフトである [EMBLmyGFF3](https://github.com/NBISweden/EMBLmyGFF3) のつくりを非常に参考にさせていただきました。
