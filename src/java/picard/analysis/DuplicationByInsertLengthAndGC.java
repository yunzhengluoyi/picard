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
 * tabulate insert size counts for both duplate and non-duplate reads
 */
@CommandLineProgramProperties(
        usage = " " +
                " ",
        usageShort = "Writes insert size counts for both duplate and non-duplate reads",
        programGroup = Metrics.class
)
public class DuplicationByInsertLengthAndGC extends SinglePassSamProgram {

    @Option(doc="If set to true, calculate mean quality over aligned reads only.")
    public boolean ALIGNED_READS_ONLY = false;

    @Option(doc="If set to true calculate mean quality over PF reads only.")
    public boolean PF_READS_ONLY = false;

    private final HistogramGenerator q  = new HistogramGenerator(false);

    /** Required main method. */
    /*public static void main(String[] args) {
        System.exit(new MeanQualityByCycle().instanceMain(args));
    }*/

    /** Required main method. */
    public static void main(String[] args) {
        System.exit(new DuplicationByInsertLengthAndGC().instanceMain(args));
    }

    private static class HistogramGenerator {
        final boolean useOriginalQualities;
        Histogram<Integer> insertHistDup = new Histogram<Integer>("GC_content", "duplicate_read_counts");
        Histogram<Integer> insertHistNotDup = new Histogram<Integer>("GC_content", "unique_read_counts");

        private HistogramGenerator(final boolean useOriginalQualities) {
            this.useOriginalQualities = useOriginalQualities;
        }

        //will hold the relevant gc information per contig
        // private int gc = -99;
        // private int gcC = -1;
        // private int referenceIndex = -1;
        // private byte [] refBases = null;
        // private java.util.ArrayList<int[]> gcCumCount = new java.util.ArrayList<int[]>();
        static int gc = -99;
        static int gcC = -1;
        static int referenceIndex = -1;
        static byte [] refBases = null;
        static java.util.ArrayList<int[]> gcCumCount = new java.util.ArrayList<int[]>();
        static int nbinsGC = 15;
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
                    if (iSz > LengthMax) {
                        nReadsOverMaxLength =+ 1; // skip read if insert length is loner than specified max
                    } else{
                        gc = GetReadGCContent(rec, ref1);
                        boolean isDup = rec.getDuplicateReadFlag();
                        //                          //
                        //  Place values in bins    //
                        //                          //
                        //System.out.println(Arrays.toString(maxOfIntvalsLen));
                        int ibinLen = (int) ((iSz - LengthMin) / binSizeLen);
                        //System.out.println(iSz + " is in bin " + maxOfIntvalsLen[ibinLen]);
                        int ibinGC = (int) Math.ceil( ((double) gc) / (binSizeGC * 100));
                        //System.out.println("gc = " + gc);
                        //System.out.println(" is in bin " + maxOfIntvalsGC[ibinGC]);
                        // add entry to result matrix
                        if (isDup) {
                            nonOptDupMtrx[ibinLen][ibinGC] += 1;
                        } else {
                            uniqueMtrx[ibinLen][ibinGC] += 1;
                        }
                        // place gc values in histogram
                        // if (isDup) {
                        //     insertHistDup.increment(Math.abs(gc));
                        // } else {
                        //     insertHistNotDup.increment(Math.abs(gc));
                        // }
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
              for (int i = 0; i < uniqueMtrx.length; i++) {
                outputWriter.write(Arrays.toString(uniqueMtrx[i]));
                outputWriter.newLine();
              }
              // write duplicate reads result matrix
              outputWriter.newLine();
              outputWriter.write("(insert size x insert GC) matrix - duplicate reads observed");
              outputWriter.newLine();
              for (int i = 0; i < uniqueMtrx.length; i++) {
                outputWriter.write(Arrays.toString(nonOptDupMtrx[i]));
                outputWriter.newLine();
              }
              outputWriter.flush();  
              outputWriter.close();  
          } catch (final IOException e) {
                System.out.println("problem writing result matrix");
          }
        }

        // public void writeMatricesToFile (String filename) {
        //     System.out.println(filename);
        //     // FileOutputStream fStream = new FileOutputStream(filename);
        //     // //ObjectOutputStream outputStream = new ObjectOutputStream(new FileOutputStream(filename));
        //     // ObjectOutputStream outputStream = new ObjectOutputStream(fStream);
        //     // outputStream.writeObject(uniqueMtrx);
        //     try {
        //         FileOutputStream fStream = new FileOutputStream(filename);
        //         ObjectOutputStream outputStream = new ObjectOutputStream(fStream);
        //         outputStream.writeObject(uniqueMtrx);
        //     } catch (final IOException e) {
        //         System.out.println("problem writting result matrix");
        //     }
        //   // for (int i = 0; i < uniqueMtrx.length; i++) {
        //   //   outputWriter.write(Arrays.toString(uniqueMtrx[i]));
        //   //   outputWriter.newLine();
            
        //   // }
        // }

    }

    @Override
    protected void setup(final SAMFileHeader header, final File samFile) {
        IOUtil.assertFileIsWritable(OUTPUT);
        q.GetMaxIntervalValues();
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
        q.printResultMtrx();
        q.writeMatricesToFile(OUTPUT.getAbsolutePath());
        System.out.println("number of reads over maxLength " + q.LengthMax + " bp: " + q.nReadsOverMaxLength);
        // final MetricsFile<?,Integer> metrics = getMetricsFile();
        // metrics.addHistogram(q.insertHistDup);
        // metrics.addHistogram(q.insertHistNotDup);
        // metrics.write(OUTPUT);
    }

}

