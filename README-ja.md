# GFF3-to-DDBJ

English version is [here](https://github.com/yamaton/gff3toddbj/blob/main/README.md).



## 概要

GFF3-to-DDBJ は、GFF3 および FASTA ファイルを、登録に必要な [DDBJ アノテーション形式](https://www.ddbj.nig.ac.jp/ddbj/file-format.html#annotation) に変換するツールです。これは `table2asn` (NCBI) や `EMBLmyGFF3` (ENA) の DDBJ 版に相当します。

出力例 (`.ann` ファイル) は [tests/golden](https://github.com/yamaton/gff3toddbj/tree/main/tests/golden) ディレクトリで確認できます。

## 精度と検証

「完璧な」GFF3 から DDBJ への変換についての厳密な定義は存在しないため、本ツールでは RefSeq の GFF3 と GenBank の対応関係をゴールドスタンダード（正解データ）としています。出力結果は以下の方法で検証しています：

1. 内部ツール `genbank-to-ddbj` を介して GenBank 由来のアノテーションと `gff3-to-ddbj` の結果を比較。
2. すべての出力を [DDBJ BioProject/BioSample/Sequence Data (MSS) Parser](https://www.ddbj.nig.ac.jp/ddbj/parser.html) に通して検証。

## インストール

### Bioconda 経由

```shell
conda create -n ddbj -c conda-forge -c bioconda gff3toddbj
conda activate ddbj
```



### PyPI 経由

```shell
conda create -n ddbj -c conda-forge -c bioconda pip samtools
conda activate ddbj
python -m pip install gff3toddbj
```



### GitHub 経由 (Nightly)

```shell
conda create -n ddbj pip
conda activate ddbj
python -m pip install 'git+https://github.com/yamaton/gff3toddbj'
```



## 使用方法

```shell
gff3-to-ddbj \
  --fasta myfile.fa \               # 必須
  --gff3 myfile.gff3 \              # 推奨 (ない場合は最低限の情報のみ生成)
  --metadata mymetadata.toml \      # 任意
  --locus_tag_prefix PREFIX_ \      # BioSample に登録がある場合は必須
  --transl_table 1 \                # デフォルト: 1 (標準)
  --output output.ann               # 任意: デフォルトは標準出力
```



### 引数の詳細

- `--locus_tag_prefix`: [BioSample で割り当てられた](https://www.ddbj.nig.ac.jp/ddbj/file-format.html#locus_tag)プレフィックス。
- `--transl_table`: 遺伝暗号表のインデックス (例: バクテリアの場合は 11)。詳細は [DDBJ 遺伝暗号表](https://www.ddbj.nig.ac.jp/ddbj/geneticcode.html) を参照してください。



## 内部処理の仕組み

GFF3-to-DDBJ は以下のパイプラインを通じてデータを処理します。

### 1. データ準備 (Data Preparation)

- **FASTA 圧縮:** 入力が標準的な Gzip 形式の場合、`bgzip` を使用して再圧縮します (例: `myfile_bgzip.fa.gz` を作成)。これによりインデックス作成が可能になり、メモリ使用量が削減されます。生成されたファイルは標準の `gzip` ツールとも互換性があります。
- **ギャップ検出:** FASTA 配列をスキャンして `N` の連なりを検出し、自動的に `assembly_gap` Feature を生成します。
- **トポロジー処理:** GFF3 に `Is_circular=true` がある場合、ツールは `TOPOLOGY` Feature を挿入し、[原点をまたぐ Feature (origin-spanning features)](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/file-formats/annotation-files/about-ncbi-gff3/#origin-spanning-features) を処理します。

### 2. Feature と Qualifier のマッピング

- **SO から INSDC への変換:** [Sequence Ontology](http://sequenceontology.org) に基づき、GFF3 の "types" を DDBJ の "Features" にマッピングします。
    - *例:* `transcript` (SO:0000673) は `misc_RNA` Feature に変換されます。
- **Qualifier のリネーム:** [リネームルール](https://github.com/yamaton/gff3toddbj/blob/main/gff3toddbj/translate_features_qualifiers.toml) に基づき、GFF3 の属性 (attributes) を DDBJ 準拠の Qualifier に変換します。
    - *例:* `ID=foobar` は `/note="ID:foobar"` になります。
- **遺伝暗号の割り当て:** ユーザー指定のインデックス (デフォルト: 1) に基づき、すべての `CDS` Feature に `/transl_table` Qualifier を自動的に追加します。

### 3. 座標処理 (Coordinate Processing)

- **結合 (Joining):** 親 (Parent) を共有する Feature は `join()` 表記を使用してマージされます。これは `CDS`、`exon`、`mat_peptide`、`V_segment`、`C_region`、`D-loop`、および `misc_feature` に適用されます。
- **RNA/Exon ロジック:** 結合された `exon` の位置情報は親となる RNA の位置情報として割り当てられ、個々の `exon` エントリは破棄されます。
    - *注意:* 直接の親が `gene` である場合、`exon` は**結合されません**。
- **部分性の処理 (Partialness):** 開始コドンまたは終止コドンが欠落している場合、`CDS` の位置情報に部分性を示す記号 (`<` または `>`) を追加します。(参照: [codon_start による翻訳開始フレームのずれ](https://www.ddbj.nig.ac.jp/ddbj/cds.html#frame))。

### 4. DDBJ 準拠ロジック (Product & Gene)

- **Product の制限:** [DDBJ の記載要領](https://www.ddbj.nig.ac.jp/ddbj/qualifiers.html#product) に準拠するため、各 `CDS` は単一の `/product` に制限されます：
    - 同一産物に複数の一般名がある場合でも、'product' に複数の名前を入力しないでください。複数の名前の区切り文字として不必要な記号を使用しないでください。2つ以上の名前を記載したい場合は、最も代表的な名前を1つ `/product` Qualifier に入力し、他は `/note` Qualifier に入力してください。
    - 名前や機能が不明な場合は、"hypothetical protein" と記載することを推奨します。
- **Gene の整合性:**
    - `/gene` Qualifier が単一の値を持つことを保証します。追加の値は `/gene_synonym` に移動されます。(参照: [Qualifier キーの定義: /gene](https://www.ddbj.nig.ac.jp/ddbj/qualifiers.html#gene))。
    - 親となる `gene` Feature から、すべての子要素 (例: `mRNA`, `CDS`) へ `/gene` および `/gene_synonym` Qualifier をコピーします。

### 5. メタデータとフィルタリング

- **メタデータの注入:** メタデータファイル から `source` 情報とグローバル Qualifier を挿入します。カスタマイズの「メタデータ設定」参照。
- **準拠フィルタリング:** [DDBJ 使用規定マトリックス](https://www.ddbj.nig.ac.jp/assets/files/pdf/ddbj/fq-j.pdf) に違反する Feature や Qualifier を削除します。
    - *注意:* このプロセスにおいて、`gene` Feature はデフォルトで破棄されます。
- **重複排除:** 処理中に生成された冗長な Qualifier 値を削除します。

### 6. 最終フォーマット

- **ソート:** 行は開始位置、[Feature の優先順位](https://github.com/yamaton/gff3toddbj/blob/1cea725cca2a8f3edb45bac45d7983e255285d5e/gff3toddbj/transforms.py#L763) (`source` と `TOPOLOGY` を最上部に配置)、および終了位置の順に並べ替えられます。

- **検証ログ:** 破棄された項目はすべて `stderr` に表示されます:

    ```
    WARNING: [Discarded] feature -------> gene (count: 49911)
    WARNING: [Discarded] (Feature, Qualifier) = (mRNA, Parent) (count: 57304)
    ```



## カスタマイズ

### メタデータ設定

GFF3/FASTA ファイルに含まれていない情報（登録者情報や共通の Qualifier など）を提供するには、TOML ファイル (例: `metadata.toml`) を使用します。

- **例:** [metadata_ddbj_example.toml](https://github.com/yamaton/gff3toddbj/blob/main/examples/metadata/metadata_ddbj_example.toml) を参照してください。
- **デフォルト:** `--metadata` が省略された場合、ツールはこの[デフォルト設定](https://github.com/yamaton/gff3toddbj/blob/main/gff3toddbj/metadata_without_COMMON.toml)を使用します。

#### 主要セクション

1. **COMMON エントリ**: `SUBMITTER`、`REFERENCE`、`COMMENT` ブロックを定義します。

2. **グローバル Qualifier (DDBJ 側の注入)**: `[COMMON.feature]` 構文を使用して、その Feature が出現するたびに特定の Qualifier を挿入するよう DDBJ システムに指示します。

    ```toml
    [COMMON.assembly_gap]
    estimated_length = "unknown"
    gap_type = "within scaffold"
    linkage_evidence = "paired-ends"
    ```

*注意: 現在、ローカル注入でサポートされているのは `[source]` と `[assembly_gap]` のみです。*



### [高度な設定] Feature と Qualifier のリネーム

GFF3 と DDBJ のフォーマットは 1:1 に対応していません。GFF3 の "types" (3列目) は DDBJ の "Features" に、GFF3 の "attributes" (9列目) は DDBJ の "Qualifiers" にマッピングされます。

`gff3-to-ddbj` は [デフォルトの変換テーブル](https://github.com/yamaton/gff3toddbj/blob/main/gff3toddbj/translate_features_qualifiers.toml) を使用してこれらを変換します。`--config_rename <FILE>` を使用してこれらのルールを上書きできます。

#### カスタマイズ例:

- **Type のリネーム:** GFF3 の type を特定の DDBJ Feature キーにマッピングします。

    ```toml
    [five_prime_UTR]
    feature_key = "5'UTR"
    ```

- **Attribute のリネーム:** GFF3 の attribute を DDBJ の Qualifier にマッピングします。`__ANY__` を使用すると、すべての Feature タイプにルールを適用できます。

    ```toml
    [__ANY__.ID]
    qualifier_key = "note"
    qualifier_value_prefix = "ID:"  # オプション
    ```

- **複雑な変換:** GFF3 の type を DDBJ の Feature/Qualifier ペアにマッピングします (例: `snRNA` を `ncRNA` に変換し、class を付与)。

    ```toml
    [snRNA]
    feature_key = "ncRNA"
    qualifier_key = "ncRNA_class"
    qualifier_value = "snRNA"
    ```

- **Attribute から Feature へのマッピング:** 特定の attribute 値を別の DDBJ Feature に変換します (例: `biotype=misc_RNA` 属性を持つ `RNA` type を `misc_RNA` Feature に変換)。

    ```toml
    [RNA.biotype.misc_RNA]
    feature_key = "misc_RNA"
    ```



### [高度な設定] Feature と Qualifier のフィルタリング

[DDBJ 使用規定マトリックス](https://www.ddbj.nig.ac.jp/assets/files/pdf/ddbj/fq-j.pdf) に準拠するため、出力は [デフォルト設定](https://github.com/yamaton/gff3toddbj/blob/main/gff3toddbj/ddbj_filter.toml) によってフィルタリングされます。この TOML ファイルで明示的に許可された Feature と Qualifier のみが最終出力に含まれます。

カスタムフィルタを使用するには、`--config_filter <FILE>` で以下の構造を持つ TOML ファイルを指定してください:

```toml
# CDS Feature にはこれらの Qualifier のみが保持されます
CDS = ["EC_number", "inference", "locus_tag", "note", "product"]
```



## トラブルシューティング

### GFF3 の検証

GFF3 ファイルを検証することをお勧めします。[GFF3 online validator](http://genometools.org/cgi-bin/gff3validator.cgi) が便利ですが、ファイルサイズは 50MB に制限されています。

### GFF3 から FASTA を分離する (必要な場合)

GFF3_to_DDBJ は、GFF3 ファイル内に `##FASTA` ディレクティブを使って FASTA 情報が含まれている場合、動作しません。付属のツール `split-fasta` は GFF3 ファイルを読み込み、GFF3 (FASTA 情報なし) と FASTA ファイルを保存します。

```shell
split-fasta path/to/myfile.gff3 --suffix "_splitted"
```

これにより、`myfile_splitted.gff3` と `myfile_splitted.fa` という2つのファイルが作成されます。



### エントリー名の正規化 (必要な場合)

DDBJ アノテーションの1列目 (= "Entry") には、`=|>" []` のような文字は使用できません。付属のプログラム `normalize-entry-names` は、そのようなエントリをリネームします。このプログラムは、例えば `ERS324955|SC|contig000013` のような ID を `ERS324955:SC:contig000013` に変換します。

```shell
normalize-entry-names myannotation_output.txt
```



このコマンドは、無効な文字が見つかった*場合*、ファイル `myannotation_output_renamed.txt` を作成します。そうでない場合、出力はありません。



## 既知の問題

### 生物学的・配列ロジック

- **トランススプライシング:** 本ツールは現在、`/trans_splicing` Qualifier を含む Feature の座標補正や `join()` 構文をサポートしていません。
- **翻訳の例外:** 開始コドンまたは終止コドンにおける `/transl_except` の座標処理はまだ実装されていません。
- **Qualifier の欠落:** `/exception` Qualifier が存在する場合でも、`/translation` Qualifier を自動生成しないため、DDBJ の検証エラーが発生する可能性があります。
- **塩基間座標:** "塩基間" の位置指定 (例: `123^124`) は現在サポートされておらず、正しく処理されない可能性があります。

### パフォーマンス

- **実行速度:** 最大限の精度を確保するため、本ツールは現在シングルプロセス構成で動作しています。大規模なゲノムデータセットでは実行時間が長くなることを想定してください。



## 謝辞

GFF3-to-DDBJ の設計は、GFF3 データを EMBL アノテーション形式に変換するための汎用ツールである [EMBLmyGFF3](https://github.com/NBISweden/EMBLmyGFF3) に触発されました。
