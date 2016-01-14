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
import picard.cmdline.CommandLineProgramTest;
import picard.annotation.RefFlatReader.RefFlatColumns;
import picard.metrics.MultilevelMetrics;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.StringUtil;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.metrics.MetricsFile;

import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Set;

/**
 * Created by hogstrom on 10/18/15.
 */
public class DuplicationByInsertLengthAndGCTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File("testdata/picard/sam/");

    public String getCommandLineProgramName() {
        return DuplicationByInsertLengthAndGC.class.getSimpleName();
    }

    @Test
    public void basic() throws Exception {
        //                                                        //
        //        Build SAM file to match GC test ref file        //
        //                                                        //
        final SAMRecordSetBuilder builder = new SAMRecordSetBuilder(true, SAMFileHeader.SortOrder.coordinate);
        builder.setRandomSeed(0);
        //final String sequence = "chr1";
        //final int sequenceIndex = builder.getHeader().getSequenceIndex(sequence);
        final int sequenceIndex = 0;
        builder.addPair("pair1", sequenceIndex, 1, 65); // expected GC = 0
        builder.addPair("pair2", sequenceIndex, 1, 165);
        builder.addPair("pair3", sequenceIndex, 1, 265);
        builder.addPair("pair4", sequenceIndex, 1, 365);
        builder.addPair("pair5", sequenceIndex, 1, 1002); //long reads should not be counted if over max

        final File samFile = File.createTempFile("tmp.DuplicationByInsertLengthAndGC.", ".sam");
        samFile.deleteOnExit();
        final SAMFileWriter samWriter = new SAMFileWriterFactory().makeSAMWriter(builder.getHeader(), false, samFile);
        for (final SAMRecord rec: builder.getRecords()) samWriter.addAlignment(rec);
        samWriter.close();

        //                                                       //
        //      Generate the metric files and run command        //
        //                                                       //
        final File outfile = File.createTempFile("test", ".GC_Length_count_matrices");
        outfile.deleteOnExit();
        final File outfileGc = File.createTempFile("test", ".insert_GC_by_dup");
        outfileGc.deleteOnExit();
        final File outfileLen = File.createTempFile("test", ".insert_Length_by_dup");
        outfileLen.deleteOnExit();

        final File ref = new File("testdata/picard/sam/GCcontentTest.fasta");
        final String[] args = new String[]{
                "INPUT=" + samFile.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "OUTPUT_GC_HIST=" + outfileGc.getAbsolutePath(),
                "OUTPUT_LEN_HIST=" + outfileLen.getAbsolutePath(),
                "REFERENCE_SEQUENCE=" + ref.getAbsolutePath(),
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);

        //                                      //
        //      Test GC Histogram values        //
        //                                      //
        final MetricsFile outputGC = new MetricsFile();
        outputGC.read(new FileReader(outfileGc));

        // Test aspects of GC histogram
        Assert.assertEquals(outputGC.getAllHistograms().size(), 1);
        final Histogram<Integer> gcHisto = outputGC.getHistogram();
        Assert.assertEquals(gcHisto.getSumOfValues(), 4.0); // check that n reads were counted (long reads excluded)
        Assert.assertEquals(gcHisto.getMin(), 0.0); // check that the insert with the lowest GC is 0

        // Test expected values of GC bin Ids
        final Set<Integer> keySetGC = gcHisto.keySet();
        final Object[] gcKeys = keySetGC.toArray();
        Assert.assertEquals(gcHisto.get(gcKeys[0]).getId(), 0); // first GC bin = 0% 
        Assert.assertEquals((int) gcHisto.get(gcKeys[0]).getValue(), 1); // should contain one value
        Assert.assertEquals(gcHisto.get(gcKeys[1]).getId(), 25); // second GC bin = 25%
        Assert.assertEquals((int) gcHisto.get(gcKeys[1]).getValue(), 1); // should contain one value
        Assert.assertEquals(gcHisto.get(gcKeys[2]).getId(), 33); // third GC bin = 33%
        Assert.assertEquals((int) gcHisto.get(gcKeys[2]).getValue(), 1); // should contain one value
        Assert.assertEquals(gcHisto.get(gcKeys[3]).getId(), 50); // third GC bin = 50%
        Assert.assertEquals((int) gcHisto.get(gcKeys[3]).getValue(), 1); // should contain one value

        //                                      //
        //     Test Length Histogram values     //
        //                                      //
        final MetricsFile outputLen = new MetricsFile();
        outputLen.read(new FileReader(outfileLen));

        // Test aspects of GC histogram
        Assert.assertEquals(outputLen.getAllHistograms().size(), 1);
        final Histogram<Integer> lenHisto = outputLen.getHistogram();
        Assert.assertEquals(lenHisto.getSumOfValues(), 4.0); // check that n reads were counted (long reads excluded)
        Assert.assertEquals( (int) lenHisto.getMin(), 100); // check that the insert with the lowest GC is 0

        // Test expected values of GC bin Ids
        final Set<Integer> lenKeysSet = lenHisto.keySet();
        final Object[] lenKeys = lenKeysSet.toArray();
        Assert.assertEquals( (int) lenHisto.get(lenKeys[0]).getId(), 100); // first GC bin = 0% 
        Assert.assertEquals( (int) lenHisto.get(lenKeys[0]).getValue(), 1); // should contain one value
        Assert.assertEquals( (int) lenHisto.get(lenKeys[1]).getId(), 200); // second GC bin = 25%
        Assert.assertEquals( (int) lenHisto.get(lenKeys[1]).getValue(), 1); // should contain one value
        Assert.assertEquals( (int) lenHisto.get(lenKeys[2]).getId(), 300); // third GC bin = 33%
        Assert.assertEquals( (int) lenHisto.get(lenKeys[2]).getValue(), 1); // should contain one value
        Assert.assertEquals( (int) lenHisto.get(lenKeys[3]).getId(), 400); // third GC bin = 50%
        Assert.assertEquals( (int) lenHisto.get(lenKeys[3]).getValue(), 1); // should contain one value

    }
}

