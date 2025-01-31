package org.go;

import java.util.ArrayDeque;
import java.util.Deque;
import java.util.HashMap;
import java.util.HashSet;

public class PathComputer
{

    private final GOEntry origin1;
    private final GOEntry origin2;
    private final HashMap<GOEntry, Integer> ancestorDistances1 = new HashMap<>();
    private final HashMap<GOEntry, Integer> ancestorDistances2 = new HashMap<>();

    public PathComputer(GOEntry term1, GOEntry term2)
    {
        this.origin1 = term1;
        this.origin2 = term2;
        exploreAncestors();
    }

    private void exploreAncestors()
    {
        ancestorDistances1.put(origin1, 0);
        ancestorDistances2.put(origin2, 0);

        performDirectedSearch(origin1, ancestorDistances1);
        performDirectedSearch(origin2, ancestorDistances2);
    }

    private void performDirectedSearch(GOEntry start, HashMap<GOEntry, Integer> distanceMap)
    {
        HashSet<GOEntry> explored = new HashSet<>();
        Deque<GOEntry> queue = new ArrayDeque<>();
        queue.add(start);

        while (!queue.isEmpty())
        {
            GOEntry current = queue.poll();
            int currentDepth = distanceMap.get(current);
            explored.add(current);

            for (GOEntry parent : current.getParents())
            {
                if (!distanceMap.containsKey(parent))
                {
                    distanceMap.put(parent, currentDepth + 1);
                    queue.add(parent);
                }
            }
        }
    }

    public int calculateShortestPath()
    {
        int minTotalDistance = Integer.MAX_VALUE;

        for (GOEntry sharedAncestor : ancestorDistances1.keySet())
        {
            if (ancestorDistances2.containsKey(sharedAncestor))
            {
                int totalDistance = ancestorDistances1.get(sharedAncestor) + ancestorDistances2.get(sharedAncestor);
                minTotalDistance = Math.min(minTotalDistance, totalDistance);
            }
        }

        return minTotalDistance == Integer.MAX_VALUE ? Integer.MAX_VALUE : minTotalDistance;
    }
}
