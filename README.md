# GFF3_to_DDBJ

## これは何？

DDBJ への登録には指定された形式のアノテーションファイルが必要です。GFF3_to_DDBJ はその名前のとおり**GFF3ファイルからアノテーションファイルを作る**プログラムです（正確には FASTA も必要です。）**FASTA 単体から最小限のアノテーションファイルを作る**ことも実は可能です。



## セットアップ

### [未実装] bioconda経由で conda環境作成

```shell
# ddbjという名前でconda環境をつくってbiocondaからパッケージをインストール
$ conda create -n ddbj -c bioconda gff3toddbj

# 環境ddbjをアクティベート
$ conda activate ddbj
```



### GitHub からダウンロードして conda 環境作成

```shell
# ファイルをダウンロード
$ wget https://github.com/yamaton/gff3_to_ddbj/archive/refs/heads/main.zip

# zipを展開
$ unzip main.zip

# ddbjという名前でconda環境をつくる
$ conda create -n ddbj

# ddbjへ environment に書かれたパッケージ群をインストール
$ conda env update -n ddbj --file environment.yaml

# 環境ddbjをアクティベート
$ conda activate ddbj
```


## アノテーションをつくる

### 0. GFF3 の正当性チェック

アノテーションファイルへの変換を始めるまえに、手持ちのファイルがGFF3形式を満たしているかチェックをかけておくのが良いと思います。[GFF3 online validator](http://genometools.org/cgi-bin/gff3validator.cgi) が便利ですが、ファイル上限が 50MB なのが玉に瑕です。



### 1. GFF3 と FASTA の分離

GFF3 ファイル中に `##FASTA` ディレクティヴをつかって FASTA が書き込まれている場合には、同梱のツールを使うなどして分割してください。`##FASTA` が無ければスキップして次へ進んでください。

```shell
split_fasta_from_gff3 \
  --gff3=path/to/myfile.gff3 \
  --suffix="_modified"
```

このばあい `myfile_modified.gff3` と `myfile_modified.fa` の2つのファイルが作られます。



### 2.（とりあえず） スクリプトの実行

まずはスクリプトを動かしてみましょう。

* `--gff3` には入力 GFF3ファイルへのパスを指定
* `--fasta` には入力 FASTAファイルへのパスを指定
* `--metadata` には４で使ったTOML ファイルを指定。（とりあえずは metadata.to
* `--locus_tag_prefix` には [BioSampleの申請時に得られたもの](https://www.ddbj.nig.ac.jp/ddbj/file-format.html#locus_tag)を指定。（とりあえずは省略してOK）
* `--transl_table` は [The Genetic Codes](https://www.ddbj.nig.ac.jp/ddbj/geneticcode-e.html) から適切な番号を選ぶ。（とりあえずは１でOK）
* `--output` には出力ファイル（＝アノテーション）のパスを指定

```shell
gff3_to_ddbj
  --gff3 myfile.gff3 \
  --fasta myfile.fa \
  --metadata metadata.toml \
  --locus_tag_prefix MYOWNPREFIX \
  --transl_table 1 \
  --output myawesome_output.ann
```

このとき出力として `myawesome_output.ann` が作成されます。



### 3. Entry 名の正規化

DDBJ のアノテーションチェックソフトによると `=|>" []` といった文字は Entry として使えないとのこと。違反文字が含まれるときには GFF3の1列目 (= "SeqID") および FASTAのヘッダをリネームする必要があります。ステップ２でその旨のエラーを見かけたら以下を実行してください。

```shell
regularize_seqids \
  --gff3=path/to/foo.gff3 \
  --fasta=path/to/bar.fasta \
  --suffix="_renamed_ids"
```

リネームの必要があるときには `foo_renamed_ids.gff3` と `bar_renamed_ids.fasta` の2つのファイルが作られます。無いときには `IDs are fine: No need to regularize them.` のメッセージが出て終了します。



### 4. メタデータの編集

DDBJアノテーションのCOMMON 項目に載せる情報や、`assembly_gap` の付属情報といったメタデータを設定するため [metadata.toml](https://raw.githubusercontent.com/yamaton/gff3_to_ddbj/main/metadata.toml) をベースに新規ファイルをつくります。テキストエディタで開いてください。COMMON を省略したいときには代わりに [metadata_without_COMMON.toml](https://raw.githubusercontent.com/yamaton/gff3_to_ddbj/main/metadata_without_COMMON.toml) から始めるのが便利です。

* COMMON に入れる[基本情報](https://www.ddbj.nig.ac.jp/ddbj/file-format.html#annotation)
* COMMON に入れる[メタ表記](https://www.ddbj.nig.ac.jp/ddbj/file-format.html#common)
* `assembly_gap` Feature に付属させる Qualifier 情報



### 5. 変換テーブルの編集

GFF3 と DDBJ アノテーションには大まかに以下のような対応があります。

1.  GFF3 の3列目 $\to$ DDBJ アノテーション 2列目 Feature
2.  GFF3 の9列目 $\to$ DDBJ アノテーション 4，5列目 Qualifier Key, Value

いっぽうでアノテーションとして許される Features と Qualifiers には制限があります。（DDBJ一覧表）

当プログラムは「よしなに」計らうよう心掛けておりますが、それでもユーザさんの手作業が要ることがあります。というのも 1, 2 における表現と DDBJの規約のあいだにはギャップがあるためです。

そのため `translate_features.toml` と `translate_qualifiers.toml` を編集することになります。



### 6.  `/product`  qualifier value の選択

Prokka などアノテーションソフトの一部は複数の値をもった `/product` を出してきます。いっぽう DDBJ アノテーションでは次のような規約があります。

> * 一般名が複数ある場合でも, 複数の名称を記載しないで下さい。また, そのために不必要な区切り記号を使用しないで下さい。一般名の複数記載を希望される場合は, 代表的な名称を /product qualifier に記載し, その他の名称を /note qualifier に記載して下さい。
>
> * 機能, 名称等が不明な蛋白質の場合は, hypothetical protein と記載することを推奨します。

このため `proteins.txt`  という各行にひとつ名称を記述したリストを準備してもらいます。ある Feature に対し`/product` が複数あるばあいには、このリストを上からたどってはじめに一致する名前をもって `/product` qualifier value とします。 リストに一致が見つからないときまたは `protein.txt`が準備されないときには `/product` qualifier value を "hypothetical protein" とします。



#### 例

protein.txt

```toml
foo
baba
```



GFF3

```
AAAA01000001.1  transdecoder    CDS 262788  262824  .   -   2   ID=cds.LOCUS000000100.1.p1;Parent=LOCUS000000100.1.p1;product=someprotein,another,baba,keke,gozilla,foo
AAAA01000001.1  transdecoder    CDS 384063  384552  .   -   0   ID=cds.LOCUS000000100.1.p1;Parent=LOCUS000000100.1.p1;product=bar,baz,abyss,nanachi
```



Annotation output (showing the corresponding part only)

|      |      |                            |             |                            |
| ---- | ---- | -------------------------- | ----------- | -------------------------- |
|      | CDS  | complement(262788..262824) | note        | ID:cds.LOCUS000000100.1.p1 |
|      |      |                            | product     | foo                        |
|      |      |                            | codon_start | 3                          |
|      | CDS  | complement(384063..384552) | note        | ID:cds.LOCUS000000100.1.p1 |
|      |      |                            | product     | hypothetical protein       |
|      |      |                            | codon_start | 1                          |





## FASTA 単体からアノテーションをつくる

GFF3+FASTA を使う前節のなかで、ステップ４のメタデータ情報の編集だけが事前に必要です。

* `--fasta` には入力 FASTAファイルへのパスを指定
* `--metadata` には４で作成したTOML ファイルを指定
* `--locus_tag_prefix` には [BioSampleの申請時に得られたもの](https://www.ddbj.nig.ac.jp/ddbj/file-format.html#locus_tag)を指定。（とりあえずは省略してOK）
* `--output` には出力ファイル（＝アノテーション）のパスを指定

```
gff3_to_ddbj
  --fasta myfile.fa \
  --metadata metadata.toml \
  --locus_tag_prefix MYOWNPREFIX \
  --output myawesome_output.ann
```
