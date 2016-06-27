package picard.sam.markduplicates;

import htsjdk.samtools.DuplicateSet;
import htsjdk.samtools.DuplicateSetIterator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.cram.encoding.NullEncoding;
import htsjdk.samtools.util.CloseableIterator;
import javafx.util.Pair;
import org.testng.annotations.Test;
import picard.PicardException;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;
import java.util.Set;
import java.util.Vector;

/**
 * Created by fleharty on 5/23/16.
 */
public class UmiAwareDuplicateSetIterator implements CloseableIterator<DuplicateSet> {
    private DuplicateSetIterator wrappedIterator;
    private Iterator<DuplicateSet> nextSetsIterator;

    public UmiAwareDuplicateSetIterator() {

    }
    public UmiAwareDuplicateSetIterator(final DuplicateSetIterator wrappedIterator) {
        this.wrappedIterator = wrappedIterator;
        nextSetsIterator = Collections.emptyIterator();
    }

    @Override
    public void close() {
        wrappedIterator.close();
    }

    @Override
    public boolean hasNext() {
        return nextSetsIterator.hasNext() || wrappedIterator.hasNext();
    }

    @Override
    public DuplicateSet next() {
        if(!nextSetsIterator.hasNext())
            process(wrappedIterator.next());

        return nextSetsIterator.next();
    }

    // Takes a duplicate set and breaks it up into possible smaller sets according to the UMI,
    // and updates nextSetsIterator to be an iterator on that set of DuplicateSets.
    private void process(final DuplicateSet set) {

        List<SAMRecord> records = set.getRecords();

        // Sort records by RX tag
        Collections.sort(records, new Comparator<SAMRecord>() {
            @Override
            public int compare(final SAMRecord lhs, final SAMRecord rhs) {
                return ((String) lhs.getAttribute("RX")).compareTo((String) rhs.getAttribute("RX"));
            }
        });

        int n = records.size();

        // Locate records that have identical UMI sequences
        int nBreaks = 1;
        List<String> uniqueObservedUMIs = new ArrayList<String>();

        uniqueObservedUMIs.add((String) records.get(0).getAttribute("RX"));
        for(int i = 1;i < n;i++) {
            // If the records differ (after sorting) we have a new duplicate set.
            if(!records.get(i).getAttribute("RX").equals(records.get(i-1).getAttribute("RX"))) {
                uniqueObservedUMIs.add((String) records.get(i).getAttribute("RX"));
                nBreaks++;
            }
            else {
            }
        }

        //for(int i = 0;i < nBreaks;i++) {
        //    System.out.println(i + " Unique UMIs " + uniqueObservedUMIs.get(i));
        //}
        int nEdges = (1+nBreaks)*nBreaks / 2;

        // Construct Adjacency List of UMIs that are close
        List<List<Integer>> adjacencyList = new ArrayList<>();
        List<Integer> groups = new ArrayList<Integer>();

        for(int i = 0;i < uniqueObservedUMIs.size();i++) {
            adjacencyList.add(i, new ArrayList<Integer>());
            groups.add(i, 0);
            for(int j = 0;j < uniqueObservedUMIs.size();j++) {
                if( getEditDistance(uniqueObservedUMIs.get(i), uniqueObservedUMIs.get(j)) <= 1) {
                    adjacencyList.get(i).add(j);
                }
            }
            //System.out.println("[ " + i + " " + adjacencyList.get(i).toString() + "]");
        }

        joinGroups(adjacencyList, groups);

//        for(int i = 0;i < groups.size();i++) {
//            System.out.println("Group " + i + " " + groups.get(i));
//        }

//        for(int i = 0;i < n;i++) {
//            System.out.println(nBreaks + " " + i + " " + records.get(i).getAttribute("RX") + " " + records.get(i).getContig() + ":" +
//                    records.get(i).getAlignmentStart() + "-"  + records.get(i).getAlignmentEnd() + " " +
//                    records.get(i).getCigarString() + " " + records.get(i).getMateAlignmentStart() + " " +
//                    records.get(i).getInferredInsertSize());
//        }

//        System.out.println(n + " " + nBreaks);

        // Construct new list of Duplicate Sets
        // Sort by group
 //       Collections.sort(groups, new Comparator<Integer>() {
 //           @Override
 //           public int compare(final Integer lhs, final Integer rhs) {
 //               return lhs.compareTo(rhs);
 //           }
 //       });

        // Figure out the number of groups
        int maxGroups = 1;
        for(int i = 0;i < groups.size();i++) {
            //System.out.println("Group " + i + " " + groups.get(i));
            if(groups.get(i) > maxGroups) {
                maxGroups = groups.get(i);
            }
        }

        // Construct DuplicateSetList
        List<DuplicateSet> duplicateSetList= new ArrayList<>();
        for(int i=0;i < maxGroups;i++) {
            DuplicateSet e = new DuplicateSet();
            duplicateSetList.add(e);
        }

        // Assign each record to a group
        for(int i = 0;i < records.size();i++) {
            String umi = (String) records.get(i).getAttribute("RX");

            // Figure out which group this belongs to
            int recordGroup = groups.get(uniqueObservedUMIs.indexOf((String) umi));
            duplicateSetList.get(recordGroup-1).add(records.get(i));
        }

        // For each group, create a duplicate set and add it to the list.
        //System.out.println("Start of Original Duplicate Set");
        for(int k = 0;k < duplicateSetList.size();k++) {
            List<SAMRecord> tmpRecords = duplicateSetList.get(k).getRecords();
            //System.out.println("Start of sub-duplicate set");
            //for (int j = 0; j < tmpRecords.size(); j++) {
            //    System.out.println("Duplicate set k = " + k + " " + tmpRecords.get(j).getAttribute("RX"));
            //}
            //System.out.println();
        }


        nextSetsIterator = duplicateSetList.iterator();

        //nextSetsIterator = Collections.singletonList(set).iterator();

    }

    List<Integer> joinGroups(List<List<Integer>> adjacencyList, List<Integer> groups) {
        int nGroups = 0;
        for(int i = 0;i < adjacencyList.size();i++) {
            // Have I visited this yet?
            if(groups.get(i) == 0) {
                // No, I haven't yet seen this
                nGroups++; // We've now seen a new group

                // Depth first search on adjacencyList, setting all the values to group
                dfs(adjacencyList, groups, i, nGroups);
            } else {
                // Yes, I have seen this before (do nothing)
            }
        }

        return groups;
    }

    void dfs(List<List<Integer>> adjacencyList, List<Integer> groups, int index, int group) {
        if(groups.get(index) == 0) {

            groups.set(index, group);

            for (int i = 0; i < adjacencyList.get(index).size(); i++) {
                dfs(adjacencyList, groups, adjacencyList.get(index).get(i), group);
            }
        } else {

        }
    }

    private int getEditDistance(String s1, String s2) {
        if(s1.length() != s2.length()) {
            throw new PicardException("Barcode lengths do not match " + s1 + " and " + s2);
        }
        int count = 0;
        for(int i = 0;i < s1.length();i++) {
            if(s1.charAt(i) != s2.charAt(i)) {
                count++;
            }
        }
        return count;
    }

    class UmiGraph {

    }

    class UmiNode {

    }

    class Edge {

    }
}

