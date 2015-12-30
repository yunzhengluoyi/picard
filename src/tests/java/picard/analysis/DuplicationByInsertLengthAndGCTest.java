/*
 * The MIT License
 *
 * Copyright (c) 2010 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package picard.analysis;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.StringUtil;
import htsjdk.samtools.metrics.MetricsFile;

import org.testng.Assert;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;
import picard.annotation.RefFlatReader.RefFlatColumns;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
/**
 * Created by hogstrom on 10/18/15.
 */
public class DuplicationByInsertLengthAndGCTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File("testdata/picard/sam/");

    public String getCommandLineProgramName() {
        return DuplicationByInsertLengthAndGC.class.getSimpleName();
    }

    @Test
    public void testGCcontent() throws IOException {
        //final File input = new File(TEST_DATA_DIR, "insert_size_metrics_test.sam");
        //final File ref = new File(TEST_DATA_DIR, "merger.fasta");
        final File ref = new File("testdata/picard/quality/chrM.reference.fasta");
        final File input = new File("testdata/picard/quality/chrMReadsDiffereingLengths.sam");
        //final File ref = new File("/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta");
        //final File input = new File(TEST_DATA_DIR, "test_chrm21_rg2.bam");
        //final File ref = new File("testdata/picard/metrics/chrMNO.reference.fasta");
        //final File input = new File("/Users/hogstrom/Documents/code/picard/testdata/picard/sam/CollectGcBiasMetrics/CollectGcBias6098159690966723109.bam");
        final File outfile = File.createTempFile("test", ".insert_GC_by_dup");
        final File pdf = File.createTempFile("test", ".pdf");
        //outfile.deleteOnExit();
        System.out.println("outfile = " +  outfile.getAbsolutePath());
        pdf.deleteOnExit();
        final String[] args = new String[]{
                "INPUT=" + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "REFERENCE_SEQUENCE=" + ref.getAbsolutePath(),
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);
        /*final MetricsFile<InsertSizeMetrics, Comparable<?>> output = new MetricsFile<InsertSizeMetrics, Comparable<?>>();
        output.read(new FileReader(outfile));

        Assert.assertEquals(output.getAllHistograms().size(), 5);

        for (final InsertSizeMetrics metrics : output.getMetrics()) {
            Assert.assertEquals(metrics.PAIR_ORIENTATION.name(), "FR");
            if (metrics.LIBRARY == null) {  // SAMPLE or ALL_READS level
                Assert.assertEquals((int) metrics.MEDIAN_INSERT_SIZE, 41);
                Assert.assertEquals(metrics.MIN_INSERT_SIZE, 36);
                Assert.assertEquals(metrics.MAX_INSERT_SIZE, 45);
                Assert.assertEquals(metrics.READ_PAIRS, 13);
                Assert.assertEquals(metrics.WIDTH_OF_10_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_20_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_30_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_40_PERCENT, 7);
                Assert.assertEquals(metrics.WIDTH_OF_50_PERCENT, 7);
                Assert.assertEquals(metrics.WIDTH_OF_60_PERCENT, 7);
                Assert.assertEquals(metrics.WIDTH_OF_70_PERCENT, 9);
                Assert.assertEquals(metrics.WIDTH_OF_80_PERCENT, 11);
                Assert.assertEquals(metrics.WIDTH_OF_90_PERCENT, 11);
                Assert.assertEquals(metrics.WIDTH_OF_99_PERCENT, 11);

            }
        }*/
    }

    @Test
    public void basic() throws Exception {
        //final String sequence = "chr1";
        final String sequence = "chrM";
        final String ignoredSequence = "chrM";

        // Create some alignments that hit the ribosomal sequence, various parts of the gene, and intergenic.
        final SAMRecordSetBuilder builder = new SAMRecordSetBuilder(true, SAMFileHeader.SortOrder.coordinate);
        // Set seed so that strandedness is consistent among runs.
        builder.setRandomSeed(0);
        final int sequenceIndex = builder.getHeader().getSequenceIndex(sequence);
        builder.addPair("pair1", sequenceIndex, 45, 475);
        builder.addPair("pair2", sequenceIndex, 90, 225);
        builder.addPair("pair3", sequenceIndex, 120, 600);
        builder.addFrag("frag1", sequenceIndex, 150, true);
        builder.addFrag("frag2", sequenceIndex, 450, true);
        builder.addFrag("frag3", sequenceIndex, 225, false);
        builder.addPair("rrnaPair", sequenceIndex, 400, 500);

        builder.addFrag("ignoredFrag", builder.getHeader().getSequenceIndex(ignoredSequence), 1, false);

        final File samFile = File.createTempFile("tmp.collectRnaSeqMetrics.", ".sam");
        //samFile.deleteOnExit();
        System.out.println("sam file = " + samFile.getAbsolutePath());

        final SAMFileWriter samWriter = new SAMFileWriterFactory().makeSAMWriter(builder.getHeader(), false, samFile);
        for (final SAMRecord rec: builder.getRecords()) samWriter.addAlignment(rec);
        samWriter.close();

        // Create an interval list with one ribosomal interval.
        final Interval rRnaInterval = new Interval(sequence, 300, 520, true, "rRNA");
        final IntervalList rRnaIntervalList = new IntervalList(builder.getHeader());
        rRnaIntervalList.add(rRnaInterval);
        final File rRnaIntervalsFile = File.createTempFile("tmp.rRna.", ".interval_list");
        rRnaIntervalsFile.deleteOnExit();
        rRnaIntervalList.write(rRnaIntervalsFile);

        // Generate the metrics.
        final File outfile = File.createTempFile("test2", ".insert_GC_by_dup");
        //metricsFile.deleteOnExit();

        /*final String[] args = new String[] {
                "INPUT=" +               samFile.getAbsolutePath(),
                "OUTPUT=" +              metricsFile.getAbsolutePath(),
                "REF_FLAT=" +            getRefFlatFile(sequence).getAbsolutePath(),
                "RIBOSOMAL_INTERVALS=" + rRnaIntervalsFile.getAbsolutePath(),
                "STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND",
                "IGNORE_SEQUENCE=" +ignoredSequence
        };*/

        final File ref = new File("testdata/picard/quality/chrM.reference.fasta");
        final String[] args = new String[]{
                "INPUT=" + samFile.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "REFERENCE_SEQUENCE=" + ref.getAbsolutePath(),
        };

        Assert.assertEquals(runPicardCommandLine(args), 0);

        /*final MetricsFile<RnaSeqMetrics, Comparable<?>> output = new MetricsFile<RnaSeqMetrics, Comparable<?>>();
        output.read(new FileReader(metricsFile));

        final RnaSeqMetrics metrics = output.getMetrics().get(0);
        Assert.assertEquals(metrics.PF_ALIGNED_BASES, 396);
        Assert.assertEquals(metrics.PF_BASES, 432);
        Assert.assertEquals(metrics.RIBOSOMAL_BASES.longValue(), 108L);
        Assert.assertEquals(metrics.CODING_BASES, 136);
        Assert.assertEquals(metrics.UTR_BASES, 51);
        Assert.assertEquals(metrics.INTRONIC_BASES, 50);
        Assert.assertEquals(metrics.INTERGENIC_BASES, 51);
        Assert.assertEquals(metrics.CORRECT_STRAND_READS, 3);
        Assert.assertEquals(metrics.INCORRECT_STRAND_READS, 4);
        Assert.assertEquals(metrics.IGNORED_READS, 1);*/
    }

    public File getRefFlatFile(String sequence) throws Exception {
        // Create a refFlat file with a single gene containing two exons, one of which is overlapped by the
        // ribosomal interval.
        final String[] refFlatFields = new String[RefFlatColumns.values().length];
        refFlatFields[RefFlatColumns.GENE_NAME.ordinal()] = "myGene";
        refFlatFields[RefFlatColumns.TRANSCRIPT_NAME.ordinal()] = "myTranscript";
        refFlatFields[RefFlatColumns.CHROMOSOME.ordinal()] = sequence;
        refFlatFields[RefFlatColumns.STRAND.ordinal()] = "+";
        refFlatFields[RefFlatColumns.TX_START.ordinal()] = "49";
        refFlatFields[RefFlatColumns.TX_END.ordinal()] = "500";
        refFlatFields[RefFlatColumns.CDS_START.ordinal()] = "74";
        refFlatFields[RefFlatColumns.CDS_END.ordinal()] = "400";
        refFlatFields[RefFlatColumns.EXON_COUNT.ordinal()] = "2";
        refFlatFields[RefFlatColumns.EXON_STARTS.ordinal()] = "49,249";
        refFlatFields[RefFlatColumns.EXON_ENDS.ordinal()] = "200,500";

        final File refFlatFile = File.createTempFile("tmp.", ".refFlat");
        //refFlatFile.deleteOnExit();
        System.out.println("refFlat = " + refFlatFile.getAbsolutePath());
        final PrintStream refFlatStream = new PrintStream(refFlatFile);
        refFlatStream.println(StringUtil.join("\t", refFlatFields));
        refFlatStream.close();

        return refFlatFile;
    }

}

