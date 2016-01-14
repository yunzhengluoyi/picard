/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
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
//package htsjdk.samtools;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.programgroups.Metrics;

import java.io.File;
import java.io.FileWriter;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.FileOutputStream;
import java.io.ObjectOutputStream;
import java.util.List;
import java.util.Arrays;

/**
 * tabulate read pair counts according insert size and insert GC for both duplate and non-duplate reads
 *
 * @author Larson Hogstrom
 */

@CommandLineProgramProperties(
        usage = " " +
                " ",
        usageShort = "Writes insert size counts and insert GC for both duplate and non-duplate read pairs. Only " +
        "the forward read of a read pair contributes to the final counts. Unmapped reads are not counted. " +
        "multiple metrics files and and (GC x insert length) count matrices. Reads with an infered insert length  " +
        "longer than 'LengthMax' will be ignored.",
        programGroup = Metrics.class
)
public class DuplicationByInsertLengthAndGC extends SinglePassSamProgram {

    @Option(doc = "File to write the output to.")
    public File OUTPUT_GC_HIST;

    @Option(doc = "File to write the output to.")
    public File OUTPUT_LEN_HIST;

    @Option(doc="If set to true, calculate counts over aligned reads only.")
    public boolean ALIGNED_READS_ONLY = false;

    @Option(doc="If set to true calculate counts over PF reads only.")
    public boolean PF_READS_ONLY = false;

    private final HistogramGenerator q  = new HistogramGenerator(false);

    /** Required main method. */
    public static void main(String[] args) {
        System.exit(new DuplicationByInsertLengthAndGC().instanceMain(args));
    }

    private static class HistogramGenerator {
        final boolean useOriginalQualities;
        // histograms for storing GC content of dups and non-dup reads
        Histogram<Integer> insertGcHistDup = new Histogram<Integer>("GC_content", "duplicate_read_counts");
        Histogram<Integer> insertGcHistNotDup = new Histogram<Integer>("GC_content", "unique_read_counts");
        // histograms for storing insert length content for dups and non-dup reads
        Histogram<Integer> insertLenHistDup = new Histogram<Integer>("Insert_length", "duplicate_reads");
        Histogram<Integer> insertLenHistNotDup = new Histogram<Integer>("Insert_length", "unique_reads");

        private HistogramGenerator(final boolean useOriginalQualities) {
            this.useOriginalQualities = useOriginalQualities;
        }

        //will hold the relevant gc information per contig
        static int gc = -99;
        static int gcC = -1;
        static int referenceIndex = -1;
        static byte [] refBases = null;
        static java.util.ArrayList<int[]> gcCumCount = new java.util.ArrayList<int[]>();
        static int nbinsGC = 30;
        static int nbinsLength = 50;
        static int LengthMin = 1;
        static int LengthMax = 1001;
        static double binSizeGC = 0;
        static int binSizeLen = 0;
        static int nReadsOverMaxLength = 0;
        //                                          //
        // Define matrices for GC and Length counts //
        //                                          //
        static int[][] uniqueMtrx = new int[nbinsLength][nbinsGC];
        static int[][] optDupMtrx = new int[nbinsLength][nbinsGC];
        static int[][] nonOptDupMtrx = new int[nbinsLength][nbinsGC];
        static double[] maxOfIntvalsGC = new double[nbinsGC]; // upper bound for bin values
        static double[] maxOfIntvalsLen = new double[nbinsLength];

        void addRecord(final SAMRecord rec, final ReferenceSequence ref1) {
            // calculate sliding cumulative GC count along reference contig
            if (!rec.getReadUnmappedFlag()) {
                boolean isForward = !rec.getReadNegativeStrandFlag();
                boolean isPaired = rec.getProperPairFlag();
                if (isForward && isPaired) {
                    int iSz = Math.abs(rec.getInferredInsertSize());
                    if (iSz >= LengthMax) {
                        nReadsOverMaxLength =+ 1; // skip read if insert length is loner than specified max
                    } else{
                        gc = GetReadGCContent(rec, ref1);
                        boolean isDup = rec.getDuplicateReadFlag();
                        //System.out.println(Arrays.toString(maxOfIntvalsLen));
                        // set bin indices
                        int ibinLen = (int) ((iSz - LengthMin) / binSizeLen);
                        int ibinGC = (int) Math.ceil( ((double) gc) / (binSizeGC * 100)) -1;
                        if (ibinGC == -1) {
                            ibinGC = 0; //if gc is content is zero set index to the first bin
                        }
                        if (ibinLen < 0 || ibinLen >= nbinsLength) {
                            System.out.println(ibinLen + ", " + iSz );
                            throw new ArrayIndexOutOfBoundsException("Length Array Index Out of Bounds");
                        }
                        if (ibinGC < 0 || ibinGC >= nbinsGC) {
                            System.out.println(ibinGC + ", " + gc );
                            throw new ArrayIndexOutOfBoundsException("GC Array Index Out of Bounds");
                        }
                        // System.out.println("gc = " + gc);
                        // System.out.println(" is in bin " + maxOfIntvalsGC[ibinGC]);
                        // System.out.println(iSz + " is in bin " + maxOfIntvalsLen[ibinLen]);
                        // add entry to result matrix
                        if (isDup) {
                            nonOptDupMtrx[ibinLen][ibinGC] += 1;
                            insertGcHistDup.increment(gc);
                            insertLenHistDup.increment(iSz);
                        } else {
                            uniqueMtrx[ibinLen][ibinGC] += 1;
                            insertGcHistNotDup.increment(gc);
                            insertLenHistNotDup.increment(iSz);
                        }
                    }
                }
            }
        }

        void printResultMtrx(){
            for (int i = 0; i < uniqueMtrx.length; i++){
                System.out.println(Arrays.toString(uniqueMtrx[i]));
            }
            System.out.println("GC intervals " + Arrays.toString(maxOfIntvalsGC));
            System.out.println("Length intervals " + Arrays.toString(maxOfIntvalsLen));
        }

        // calculate sliding cumulative GC count along reference contig
        void GetContigGCCumulativeSum (final SAMRecord rec, final ReferenceSequence ref1) {
            if (referenceIndex != rec.getReferenceIndex() || gc == -99) {
                System.out.println("calculating ref index i=" + rec.getReferenceIndex());
                refBases = ref1.getBases();
                StringUtil.toUpperCase(refBases);
                final int refLength = refBases.length;
                int gcCount = 0;
                int nCount = 0;
                gcCumCount.clear();
                int[] gcCumList = new int[refLength];
                for (int i = 0; i < refLength; ++i) {
                    final byte base = refBases[i];
                    if (SequenceUtil.basesEqual(base, (byte) 'G') || SequenceUtil.basesEqual(base, (byte) 'C')) {
                        ++gcCount;
                    } else if (SequenceUtil.basesEqual(base, (byte) 'N')) {
                        ++nCount;
                    }
                    gcCumList[i] = gcCount;
                }
                gcCumCount.add(gcCumList);
                referenceIndex = rec.getReferenceIndex();
            }
        }

        int GetReadGCContent (final SAMRecord rec, final ReferenceSequence ref1) {
                GetContigGCCumulativeSum(rec, ref1);
                int iSz = Math.abs(rec.getInferredInsertSize());
                int iStart = rec.getStart() - 1;
                //System.out.println("iStart = " + iStart + ", iSize = " + iSz);
                final int[] gcL = gcCumCount.get(0);
                boolean outBounds = (iStart < 0) || (iSz < 0) || ((iStart  +iSz) >= gcL.length);
                final int szMax = 5000;
                if (outBounds || (iSz > szMax) || (gcL.length <=0)) {
                    gcC = -1;
                } else {
                    gcC = (gcL[iStart + iSz - 1] - gcL[iStart]);
                    //System.out.println("cumsum Start = " + gcL[iStart + iSz] + ", stop = " + gcL[iStart]);
                }
                gc = -99;
                if (iSz != 0) {
                    gc = (gcC * 100) / (iSz);
                }
                return gc;
        }

        // Obtain the max values for bins - GC and insert size 
        public void GetMaxIntervalValues () {
            //                                      //
            //  Set intervals for insert length     //
            //                                      //
            binSizeGC = 1.0/ (double)(nbinsGC);
            // create array with max value for each interval
            double stepGC = 0.0;
            for (int i = 0; i < nbinsGC; i++){
                stepGC = stepGC + binSizeGC;
                maxOfIntvalsGC[i] = stepGC; 
            }
            //                                      //
            //  Set intervals for insert length     //
            //                                      //
            if (((LengthMax - LengthMin) % nbinsLength) != 0) {
                System.out.println("bin interval for insert length must result in integer values");
                System.exit(1);
            }
            else {
                binSizeLen = (LengthMax - LengthMin) / nbinsLength; 
            }
            // create array with max value for each interval
            double stepLen = LengthMin;
            for (int i = 0; i < nbinsLength; i++){
                stepLen = stepLen + binSizeLen;
                maxOfIntvalsLen[i] = stepLen; // upper bound for bin values
            }
        }

        //public void writeMatricesToFile (String filename) throws IOException {
        public void writeMatricesToFile (String filename) {    
          try {
              BufferedWriter outputWriter = null;
              outputWriter = new BufferedWriter(new FileWriter(filename));
              // write unique reads result matrix
              outputWriter.write("(insert size x insert GC) matrix - Unique reads observed");
              outputWriter.newLine();
              outputWriter.write("-,");
              outputWriter.write(Arrays.toString(maxOfIntvalsGC)
                        .replace("[","")
                        .replace("]",""));
              outputWriter.newLine();
              for (int i = 0; i < uniqueMtrx.length; i++) {
                outputWriter.write(maxOfIntvalsLen[i] + ", ");
                outputWriter.write(Arrays.toString(uniqueMtrx[i])
                        .replace("[","")
                        .replace("]",""));
                outputWriter.newLine();
              }
              // write duplicate reads result matrix
              outputWriter.newLine();
              outputWriter.write("(insert size x insert GC) matrix - duplicate reads observed");
              outputWriter.newLine();
              outputWriter.write("-,");
              outputWriter.write(Arrays.toString(maxOfIntvalsGC)
                        .replace("[","")
                        .replace("]",""));
              outputWriter.newLine();
              for (int i = 0; i < uniqueMtrx.length; i++) {
                outputWriter.write(maxOfIntvalsLen[i] + ", ");
                outputWriter.write(Arrays.toString(nonOptDupMtrx[i])
                        .replace("[","")
                        .replace("]",""));
                outputWriter.newLine();
              }
              outputWriter.flush();  
              outputWriter.close();  
          } catch (final IOException e) {
                System.out.println("problem writing result matrix");
          }
        }
    }

    @Override
    protected void setup(final SAMFileHeader header, final File samFile) {
        IOUtil.assertFileIsWritable(OUTPUT);
        q.GetMaxIntervalValues();
        //q.printResultMtrx();
    }

    @Override
    protected void acceptRead(final SAMRecord rec, final ReferenceSequence ref) {
        // Skip unwanted records
        if (PF_READS_ONLY && rec.getReadFailsVendorQualityCheckFlag()) return;
        if (ALIGNED_READS_ONLY && rec.getReadUnmappedFlag()) return;
        if (rec.isSecondaryOrSupplementary()) return;
        q.addRecord(rec,ref);
    }

    @Override
    protected void finish() {
        // Generate a "Histogram" of insert size length
        //q.printResultMtrx();
        q.writeMatricesToFile(OUTPUT.getAbsolutePath());
        //System.out.println("number of reads over maxLength " + q.LengthMax + " bp: " + q.nReadsOverMaxLength);
        final MetricsFile<?,Integer> metrics1 = getMetricsFile();
        metrics1.addHistogram(q.insertGcHistDup);
        metrics1.addHistogram(q.insertGcHistNotDup);
        metrics1.write(OUTPUT_GC_HIST);
        final MetricsFile<?,Integer> metrics2 = getMetricsFile();
        metrics2.addHistogram(q.insertLenHistDup);
        metrics2.addHistogram(q.insertLenHistNotDup); 
        metrics2.write(OUTPUT_LEN_HIST);       
    }

}

