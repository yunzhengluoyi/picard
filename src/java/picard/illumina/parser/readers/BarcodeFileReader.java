package picard.illumina.parser.readers;

import htsjdk.samtools.util.CloseableIterator;
import picard.util.BasicInputParser;

import java.io.File;

/**
 * Reads a single barcode file line by line and returns the barcode if there was a match or NULL otherwise.
 *
 * Barcode.txt file Format (consists of tab delimited columns, 1 record per row)
 * sequence_read    Matched(Y/N)    BarcodeSequenceMatched
 *
 * sequence read          - the actual bases at barcode position
 * Matched(y/n)           - Y or N indicating if there was a barcode match
 * BarcodeSequenceMatched - matched barcode sequence (empty if read did not match one of the barcodes).
 */
public class BarcodeFileReader extends AbstractIlluminaFileReader implements CloseableIterator<String> {
    private static final int Y_OR_N_COLUMN = 1;
    private static final int BARCODE_COLUMN = 2;
    private final BasicInputParser textIterator;
    private long numClusters = 0;

    public BarcodeFileReader(final File barcodeFile) {
        this(barcodeFile, false);
    }

    public BarcodeFileReader(final File barcodeFile, boolean count) {
        this.textIterator = new BasicInputParser(false, barcodeFile);
        if (count) {
            BasicInputParser counter = new BasicInputParser(false, barcodeFile);
            while (counter.hasNext()) {
                counter.next();
                numClusters++;
            }
            counter.close();
        }
    }

    @Override
    public String next() {
        final String [] fields = textIterator.next();
        final String barcode;
        if (fields[Y_OR_N_COLUMN].equals("Y")) {
            barcode = fields[BARCODE_COLUMN];
        } else {
            barcode = null;
        }

        return barcode;
    }

    @Override
    public boolean hasNext() {
        return textIterator.hasNext();
    }

    public void remove() {
        throw new UnsupportedOperationException("Remove is not supported by " + BarcodeFileReader.class.getName());
    }

    public void close() {
        textIterator.close();
    }

    @Override
    public Long getNumClusters() {
        return numClusters;
    }
}
