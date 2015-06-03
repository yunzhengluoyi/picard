package picard.illumina.parser;

import picard.illumina.parser.fakers.FileFaker;
import picard.illumina.parser.readers.AbstractIlluminaFileReader;
import picard.illumina.parser.readers.IlluminaFileReaderFactory;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

public class PerTileFileUtil extends ParameterizedFileUtil {
    private final IlluminaFileMap fileMap;

    public PerTileFileUtil(final String extension, final File base,
                           final FileFaker faker, final int lane) {
        super(true, extension, base, faker, lane);
        this.fileMap = getTiledFiles(base, matchPattern);
        if (fileMap.size() > 0) {
            this.tiles = Collections.unmodifiableList(new ArrayList<Integer>(this.fileMap.keySet()));
        } else {
            this.tiles = new ArrayList<Integer>();
        }
    }

    @Override
    public boolean filesAvailable() {
        return !fileMap.isEmpty();
    }

    public IlluminaFileMap getFiles() {
        return fileMap;
    }

    public IlluminaFileMap getFiles(final List<Integer> tiles) {
        return fileMap.keep(tiles);
    }

    @Override
    public List<String> verify(final Map<Integer, Long> expectedTiles, final int[] expectedCycles, IlluminaFileUtil.SupportedIlluminaFormat format,
                               boolean fullTest) {
        final List<String> failures = new LinkedList<String>();

        if (!base.exists()) {
            failures.add("Base directory(" + base.getAbsolutePath() + ") does not exist!");
        } else {
            if (!tiles.containsAll(expectedTiles.keySet())) {
                final List<Integer> missing = new ArrayList<Integer>(expectedTiles.keySet());
                    missing.removeAll(tiles);
                    failures.add("Missing tile " + missing + " for file type " + extension + ".");
            } else if (fullTest) {
                for (Integer tile : fileMap.keySet()) {
                    File file = fileMap.get(tile);
                    AbstractIlluminaFileReader reader = IlluminaFileReaderFactory.makeReader(format, file);
                    long numClusters = reader.getNumClusters();
                    if (numClusters != expectedTiles.get(tile)) {
                        failures.add("Expected " + expectedTiles.get(tile) + " for file " + file.getAbsolutePath() + " but only found " + numClusters);
                    }
                }
                }

        }
        return failures;
    }

    @Override
    public List<String> fakeFiles(final List<Integer> expectedTiles, final int[] cycles,
                                  final IlluminaFileUtil.SupportedIlluminaFormat format) {
        final List<String> failures = new LinkedList<String>();
        if (!base.exists()) {
            failures.add("Base directory(" + base.getAbsolutePath() + ") does not exist!");
        } else {
            for (final Integer tile : expectedTiles) {
                if (!tiles.contains(tile) || fileMap.get(tile).length() == 0) {
                    //create a new file of this type
                    try {
                        faker.fakeFile(base, tile, lane, extension);
                    } catch (final IOException e) {
                        failures.add(String.format("Could not create fake file %s: %s", fileMap.get(tile),
                                e.getMessage()));
                    }

                }
            }
        }
        return failures;
    }

}
