# GFF3-to-DDBJ

[TOC]

## これは何？

DDBJ への登録には指定された形式のアノテーションファイルが必要です。GFF3-to-DDBJ は **FASTA と GFF3ファイルからアノテーションファイルを作る**プログラムです。**FASTA 単体から最小限のアノテーションファイルを作る**ことも可能です。



## セットアップ

#### [未登録] biocondaからconda環境にインストールする場合

```shell
# ddbjという名前でconda環境をつくってbiocondaからパッケージをインストール
## 現在登録処理中のためごちゃついてます 2021-09-14
$ conda create -n ddbj -c bioconda -c conda-forge -c https://168588-42372094-gh.circle-artifacts.com/0/tmp/artifacts/packages gff3toddbj

# 環境ddbjをアクティベート
$ conda activate ddbj
```



#### GitHubソースコードからconda環境にインストールする場合

```shell
# ファイルをダウンロード
$ wget https://github.com/yamaton/gff3_to_ddbj/archive/refs/heads/main.zip

# zipを展開、リネーム
$ unzip main.zip && mv gff3toddbj-main gff3toddbj && cd gff3toddbj

# ddbjという名前でconda環境をつくる
$ conda create -n ddbj

# 環境ddbjをアクティベート
$ conda activate ddbj

# ddbjへ 依存パッケージ (bioconda, bcbio-gff, toml) をインストール
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
  --metadata mymetadata.toml \       # この行を削除するとそれなりの出力に
  --locus_tag_prefix MYOWNPREFIX \   # この行を削除すると LOCUSTAGPREFIX_ に設定
  --transl_table 1 \                 # この行を削除すると 1 と設定
  --output myawesome_output.ann      # この行を削除すると標準出力に
```



## 設定いろいろ

### 設定ファイル `config.toml` の編集

より良い出力のため設定ファイルをコピペ＆編集することをお勧めします。DDBJアノテーションのCOMMON 項目に載せる情報や、`assembly_gap` の付属情報といったを設定するため [config.toml](https://raw.githubusercontent.com/yamaton/gff3toddbj/main/gff3toddbj/config.toml) をベースに新規ファイルをつくります。テキストエディタで開いてください。COMMON を自分で追加するばあいには代わりに [config_without_COMMON.toml](https://raw.githubusercontent.com/yamaton/gff3toddbj/main/gff3toddbj/config_without_COMMON.toml) から始めるのが便利です。

* COMMON に入れる[基本情報](https://www.ddbj.nig.ac.jp/ddbj/file-format.html#annotation)

* COMMON に入れる[メタ表記](https://www.ddbj.nig.ac.jp/ddbj/file-format.html#common)

  * COMMON の下に Feature を入れておくことで、DDBJ でつくられるフラットファイルに自動的に一律に情報が挿入される機能があるそうです。たとえば以下のようにしておくと `assembly_gap` ごとに同じ Qualifier値が挿入されることになります。

    ```toml
    [COMMON.assembly_gap]
    estimated_length = "unknown"
    gap_type = "within scaffold"
    linkage_evidence = "paired-ends"
    ```

* Feature ごとに挿入する Qualifier 情報

  * Feature ごとのQualifier値の一律挿入を **gff3-to-ddbj がアノテーションファイルに対して**行うのがこの設定値です。「COMMONに入れるメタ情報」と実質的に同じことですが、DDBJで行われるメタ表記の挙動が未確認なのでこの機能を付けています。使い方は`[COMMON.assembly_gap]` を `[assembly_gap]` に置き換えるだけです。

    ```toml
    [assembly_gap]
    estimated_length = "unknown"
    gap_type = "within scaffold"
    linkage_evidence = "paired-ends"
    ```



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
  --config config.toml \
  --locus_tag_prefix MYOWNPREFIX_ \
  --transl_table 1 \
  --translate_features translate_features.toml \
  --translate_qualifiers  translate_qualifiers.toml \
  --output myawesome_output.ann
```



## トラブルシューティング

### GFF3 の正当性チェック

アノテーションファイルへの変換を始めるまえに、手持ちのファイルがGFF3形式を満たしているかチェックをかけておくのが良いと思います。[GFF3 online validator](http://genometools.org/cgi-bin/gff3validator.cgi) が便利ですが、ファイル上限が 50MB なのが玉に瑕です。



### GFF3 と FASTA の分離（必要に応じて）

GFF3 ファイル中に `##FASTA` ディレクティヴをつかって FASTA が書き込まれている場合には、同梱のツールを使うなどして分割してください。`##FASTA` が無ければスキップして次へ進んでください。

```shell
split-fasta path/to/myfile.gff3 --suffix "_modified"
```

このばあい `myfile_modified.gff3` と `myfile_modified.fa` の2つのファイルが作られます。



### Entry 名の正規化（必要に応じて）

DDBJ のアノテーションチェックソフトによると `=|>" []` といった文字は Entry として使えないとのことです。違反文字が含まれるときには GFF3の1列目 (= "SeqID") および FASTAのヘッダをリネームする必要があります。ステップ２でその旨のエラーを見かけたら以下を実行してください。

```shell
rename-ids \
  --gff3=path/to/foo.gff3 \     # <<必須>>
  --fasta=path/to/bar.fasta \   # <<必須>>
  --suffix="_renamed"       # この行を省略するとデフォルト値に
```

リネームの必要があるときには `foo_renamed.gff3` と `bar_renamed.fasta` の2つのファイルが作られます。無いときには `IDs are fine: No need to regularize them.` のメッセージが出るだけで終了します。



## 当プログラムが行うこと

表示形式を変えるほかに以下のようなことをしています。

* 変換テーブルに基づいた Features / Qualifiers のリネーム

* assembly_gap の検索

* /transl_table の CDS への追加

* メタデータ中の source 情報を各エントリへの追加

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

