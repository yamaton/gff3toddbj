# `gff3-to-ddbj` の評価



## RefSeqデータを使っての比較

![](../images/evaluation.png)

RefSeqが GFF3+FASTA および GenBank 形式を公開しているので、これを利用して `gff3-to-ddbj` の評価をしています。`gff3-to-ddbj` はGFF3とFASTAから、`genbank-to-ddbj`はGenBank形式からそれぞれDDBJアノテーションを作成したうえで、2ファイルを比べています。ここで `genbank-to-ddbj` はGenbank形式の読み取りとDDBJ形式への書き込みに限定したつくりのため、評価に影響はないものと仮定しています。

評価は２つのDDBJ出力に対して現在３点で比べることで行います。どの場合も[多重集合](https://ja.wikipedia.org/wiki/多重集合) (multiset) を作ったうえで、交わりと全体の要素数比から左側・右側それぞれについて**不一致**具合を求めています。ベン図でいう以下のような分子分母の計算です。

![](../images/fraction.png)



#### 1. 位置情報（開始・終了点補正を無視）

アノテーションファイルの位置情報を含む各行から（エントリ名、Feature名、Location補正なし）を要素とする多重集合を作ります。Locationでは開始点・終了点が不明であることを示す不等号は無視します。たとえば `join(12..78,134..202)` と `join(<12..78,134..202)`は同じものとみなします。



#### 2. 位置情報（補正等すべて含む）

１と同様に（エントリ名、Feature名、Location情報）を要素とする多重集合をつくります。ここでlocation情報には不等号等も含めます。つまり`join(12..78,134..202)` と `join(<12..78,134..202)`は違うものとして扱います。



### 3. Feature-Qualifier 情報

アノテーションファイルの位置情報を含む各行から（エントリ名、Feature名、Qualifierキー名）を要素とする多重集合を作ります。





## 評価の手順

1. RefSeq から評価用データをダウンロード (全配列 .fasta、アノテーション .gff、配列＋アノテーション .gbff の3種)

2. それぞれの tarファイルを展開して .gz ファイルを取り出しておく

3. GenBank から DDBJ に変換

   ```shell
   genbank-to-ddbj --gbk ./GCF_017312705.1_Mj_TUMSAT_v1.0_genomic.gbff.gz > refseq_genbank.ann 2> refseq_genbank.err
   ```

4. FASTA + GFF から DDBJ に変換

   ```shell
   gff3-to-ddbj --fasta ./foo.fna.gz --gff3 ./foo.gff.gz --output refseq_gff.ann 2> refseq_gff.err
   ```

5. ２つの DDBJ 形式を比較

   ```shell
   compare-ddbj refseq_gff.ann refseq_genbank.ann
   ```

   `compare-ddbj` コマンドは、それぞれの多重集合データをタブ区切りテキストとして書き出します。



## RefSeqのマウスの参照ゲノムでの結果

### サマリー
[RefSeqの代表的なゲノム](https://www.ncbi.nlm.nih.gov/genome/annotation_euk/all/?utm_source=blog&utm_medium=referrer&utm_campaign=gdv&utm_term=intron&utm_content=20210202link1)のひとつである、マウスゲノム(GRCm39) についての評価サマリーです。

```shell
$ compare-ddbj GCF_000001635.27_GRCm39_genomic.ann refseq_GCF_000001635.27_GRCm39_genomic.ann
Stat w/o location correction:
    Left  mismatching: 1248 / 238958    (0.52 %)
    Right mismatching: 575 / 238285     (0.24 %)
Stat with location correction:
    Left  mismatching: 1988 / 238958    (0.83 %)
    Right mismatching: 1315 / 238285    (0.55 %)
Stat of feature-qualifier pairs:
    Left  mismatching: 104990 / 829632  (12.66 %)
    Right mismatching: 292303 / 1016945 (28.74 %)
```


### 左側 Feature-Qualifier の差異
```shell
$ cat quals_left-only.txt | awk '{print $2, $3}' | sort | uniq -c | sort -rhb
  93152 CDS transl_table
   8302 exon pseudogene
    666 exon note
    666 exon gene
    666 CDS product
    530 misc_RNA pseudogene
    268 assembly_gap linkage_evidence
    176 CDS pseudogene
    153 V_segment pseudogene
    134 assembly_gap gap_type
    134 assembly_gap estimated_length
     63 source note
     21 tRNA note
     12 J_segment pseudogene
     10 regulatory note
     10 D_segment pseudogene
      8 protein_bind note
      7 CDS note
      6 CDS exception
      4 ncRNA_lncRNA pseudogene
      1 rRNA note
      1 C_region pseudogene
```


### 右側 Feature-Qualifier の差異分

```
$ cat quals_right-only.txt | awk '{print $2, $3}' | sort | uniq -c | sort -rhb
  92499 CDS translation
  86831 CDS gene_synonym
  86264 mRNA gene_synonym
   9456 misc_RNA gene_synonym
   8126 exon pseudo
   2409 exon gene_synonym
   2112 ncRNA_miRNA gene_synonym
   1353 ncRNA_lncRNA gene_synonym
   1228 precursor_RNA gene_synonym
    530 misc_RNA pseudo
    431 V_segment gene_synonym
    176 CDS pseudo
    153 V_segment pseudo
    118 regulatory regulatory_class
    118 regulatory experiment
    110 ncRNA_snoRNA gene_synonym
    100 regulatory note
     90 J_segment gene_synonym
     31 rRNA gene_synonym
     24 D_segment gene_synonym
     21 C_region gene_synonym
     13 ncRNA_snRNA gene_synonym
     13 CDS function
     13 CDS EC_number
     12 J_segment pseudo
     10 D_segment pseudo
      9 ncRNA_guide_RNA gene_synonym
      9 assembly_gap gap_type
      9 assembly_gap estimated_length
      6 ncRNA_antisense_RNA gene_synonym
      6 CDS ribosomal_slippage
      4 ncRNA_RNase_P_RNA gene_synonym
      4 ncRNA_lncRNA pseudo
      3 regulatory gene_synonym
      2 ncRNA_scRNA gene_synonym
      2 misc_feature note
      2 misc_feature gene_synonym
      2 misc_feature gene
      1 source organelle
      1 ncRNA_telomerase_RNA gene_synonym
      1 ncRNA_RNase_MRP_RNA gene_synonym
      1 C_region pseudo
```

## 差異の原因

### trivial な点
* GenBank形式はフラットファイルとされる最終形式に対して、DDBJアノテーションファイルはフラットファイルに変換される以前の形式。そのため
    * `/translation` は DDBJアノテーションには含まれない
* DDBJ では `assembly_gap` をきちんと付ける傾向がありますが、RefSeqはそうでもないようです。
* DDBJ では `/pseudo` を新規登録では使わない方針を示していますが、RefSeqはそうでもないようです。
* DDBJ では `transl_table` を各CDSに付ける傾向にありますが、RefSeqはそうでもないようです。


### 改善の可能性として

* RefSeq は 各種RNA と CDS Feature が locationにおいて完全一致するばあい GenBank出力にて RNAを削除することがあります。
    * 同様の操作がDDBJ登録で求められるのか確認中。
* NCBI には < 10bp 長のイントロンを載せない方針があります。結果として、`exon` 間の隙間が小さいときには `join()`記法ではなく、ひとつの`exon` としてGenBank版に記載されています。
    * 同様の操作がDDBJ登録で求められるのか確認中。

* `gff3-to-ddbj` は CDS にのみ、start codon, stop codon に応じた `<` や `>` による開始・終了点の補正を行っています。いっぽう RefSeq の.gb形式は各種RNAなど、CDS の上位階層のFeatureにもCDSに対応した位置の補正が行われています。
    * 同様の操作がDDBJ登録で求められるのか確認中。

* RefSeq アノテーションプログラム (gnomon) は入力データに無い塩基対を補足することがあり、そのばあい以下のような `/note` が付けられます。これにより開始・終了点補正がさまざまなところに入れられて、locationが異なるという結果になります。
    * これは変換ツールとしての `gff3-to-ddbj`の守備範囲外とします

  ```
         /note="The sequence of the model RefSeq protein was
          modified relative to this genomic sequence to represent
          the inferred CDS: added 170 bases not found in genome
          assembly; Derived by automated computational analysis
          using gene prediction method: Gnomon."
  ```

* `gff3-to-ddbj` は start codon, stop codon を Biopython の `Bio.Data.CodonTable.unambiguous_dna_by_id[transl_table]` に基づいて一律に判定している反面、RefSeq は状況に応じた start/stop codon 判定をしているようです。たとえば Genetic Code が 1 のときに CTG や TTG  が start codon として扱われない場合が確認されています。
    * 作者が生物知識に疎いため深掘り出来ていません。

* RefSeq は `gene` feature における `/gene` や `/gene_synonym` を下の階層の qualifiers にもコピーする。
    * 導入検討中

