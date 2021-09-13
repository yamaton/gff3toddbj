# GFF3_to_DDBJ

[TOC]

## これは何？

DDBJ への登録には指定された形式のアノテーションファイルが必要です。GFF3_to_DDBJ はその名前のとおり**GFF3ファイルからアノテーションファイルを作る**プログラムです（正確には FASTA も必要です。）**FASTA 単体から最小限のアノテーションファイルを作る**ことも実は可能です。



## セットアップ

#### bioconda経由で conda環境をつくる場合

```shell
# ddbjという名前でconda環境をつくってbiocondaからパッケージをインストール
$ conda create -n ddbj -c bioconda gff3toddbj

# 環境ddbjをアクティベート
$ conda activate ddbj
```



#### GitHub からダウンロードして conda環境で動かす場合

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



## GFF3 と FASTA から DDBJ アノテーションをつくる

### 0. GFF3 の正当性チェック

アノテーションファイルへの変換を始めるまえに、手持ちのファイルがGFF3形式を満たしているかチェックをかけておくのが良いと思います。[GFF3 online validator](http://genometools.org/cgi-bin/gff3validator.cgi) が便利ですが、ファイル上限が 50MB なのが玉に瑕です。



### 1. GFF3 と FASTA の分離

GFF3 ファイル中に `##FASTA` ディレクティヴをつかって FASTA が書き込まれている場合には、同梱のツールを使うなどして分割してください。`##FASTA` が無ければスキップして次へ進んでください。

```shell
tools/split_fasta_from_gff3 \
  --gff3=path/to/myfile.gff3 \
  --suffix="_modified"
```

このばあい `myfile_modified.gff3` と `myfile_modified.fa` の2つのファイルが作られます。



### 2.（とりあえず） スクリプトの実行

まずはスクリプトを動かしてみましょう。

* `--gff3 <FILE>` には入力 GFF3ファイルへのパスを指定
* `--fasta <FILE>` には入力 FASTAファイルへのパスを指定
* `--metadata <FILE>` には４でコピー＆編集するTOML ファイルを指定。（とりあえずはテンプレートを使ってみます。）
* `--locus_tag_prefix <STRING>` には [BioSampleの申請時に得られたもの](https://www.ddbj.nig.ac.jp/ddbj/file-format.html#locus_tag)を指定。（とりあえずは省略）
* `--transl_table <INT>` は [The Genetic Codes](https://www.ddbj.nig.ac.jp/ddbj/geneticcode-e.html) から適切な数字を選ぶ。（とりあえずは１）
* `--output <FILE>` には出力ファイル（＝アノテーション）のパスを指定

```shell
# gff3_to_ddbj のような実行ファイルを用意する予定です…
python gff3toddbj/main.py
  --gff3 myfile.gff3 \
  --fasta myfile.fa \
  --metadata mymetadata.toml \
  --locus_tag_prefix MYOWNPREFIX \
  --transl_table 1 \
  --output myawesome_output.ann
```

このとき出力として `myawesome_output.ann` が作成されます。



### 3. Entry 名の正規化

DDBJ のアノテーションチェックソフトによると `=|>" []` といった文字は Entry として使えないとのこと。違反文字が含まれるときには GFF3の1列目 (= "SeqID") および FASTAのヘッダをリネームする必要があります。ステップ２でその旨のエラーを見かけたら以下を実行してください。

```shell
tools/regularize_seqids \
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





## FASTA単体 から DDBJ アノテーションをつくる

GFF3+FASTA を使う前節のなかで、ステップ４のメタデータ情報の編集だけが事前に必要です。

* `--fasta` には入力 FASTAファイルへのパスを指定
* `--metadata` には４で作成したTOML ファイルを指定
* `--locus_tag_prefix` には [BioSampleの申請時に得られたもの](https://www.ddbj.nig.ac.jp/ddbj/file-format.html#locus_tag)を指定。（とりあえずは省略してOK）
* `--output` には出力ファイル（＝アノテーション）のパスを指定

```
python gff3toddbj/main.py
  --fasta myfile.fa \
  --metadata metadata.toml \
  --locus_tag_prefix MYOWNPREFIX \
  --output myawesome_output.ann
```



## 当プログラムが行っていること

表示形式を変えるほかに以下のようなことをしています。

* 変換テーブルに基づいた Features / Qualifiers のリネーム

* assembly_gap の検索

* /transl_table の CDS への追加

* メタデータ中の source 情報を各エントリへの追加

* Qualifier 値に使われる文字の正規化

* 同じ親をもつ CDS を `join` 記法で結合

* 入力GFF3 ファイル中の mRNA と exon を `join`記法で結合し mRNA Feature とする

* Start codon まわり整合性チェック（ いまのところ /codon_start が 1でないときのみ）

* CDS 下の /product が値をひとつだけ持つよう変更。複数値の残りは /note へ

  * [DDBJ の/product 詳細](https://www.ddbj.nig.ac.jp/ddbj/qualifiers.html#product)を参照ください

    > * 一般名が複数ある場合でも, 複数の名称を記載しないで下さい。また, そのために不必要な区切り記号を使用しないで下さい。一般名の複数記載を希望される場合は, 代表的な名称を /product qualifier に記載し, その他の名称を /note qualifier に記載して下さい。
    > * 機能, 名称等が不明な蛋白質の場合は, hypothetical protein と記載することを推奨します。

* Qualifier値に重複があるときに冗長分を削除

* アノテーション行の並び替え

* [DDBJ の Feature-Qualifier 一覧表](https://www.ddbj.nig.ac.jp/assets/files/pdf/ddbj/fq-j.pdf)に基づいた出力情報のフィルタリング

