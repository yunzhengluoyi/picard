package picard.illumina.parser.readers;

import picard.PicardException;
import picard.illumina.parser.IlluminaFileUtil;

import java.io.File;

public class IlluminaFileReaderFactory {
    public static AbstractIlluminaFileReader makeReader(IlluminaFileUtil.SupportedIlluminaFormat format, File file) {
        AbstractIlluminaFileReader reader;
        switch (format) {
            case Barcode:
                reader = new BarcodeFileReader(file, true);
                break;

            case Bcl:
                reader = new BclReader(file, new BclQualityEvaluationStrategy(BclQualityEvaluationStrategy.ILLUMINA_ALLEGED_MINIMUM_QUALITY), false);
                break;

            case Filter:
                reader = new FilterFileReader(file);
                break;

            case Locs:
                reader = new LocsFileReader(file);
                break;

            case Clocs:
                reader = new ClocsFileReader(file, true);
                break;

            case Pos:
                reader = new PosFileReader(file, true);
                break;

            default:
                throw new PicardException("Unrecognized data type(" + format + ") found by IlluminaDataProviderFactory!");
        }
        return reader;
    }
}
