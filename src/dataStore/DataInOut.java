package dataStore;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import static utilities.Utilities.*;
import com.bbn.openmap.omGraphics.OMGraphic;
import com.bbn.openmap.dataAccess.shape.EsriPolyline;
import com.bbn.openmap.dataAccess.shape.EsriPolylineList;
import com.bbn.openmap.dataAccess.shape.EsriShapeExport;
import com.bbn.openmap.dataAccess.shape.DbfTableModel;
import com.bbn.openmap.dataAccess.shape.EsriPoint;
import com.bbn.openmap.dataAccess.shape.EsriPointList;
import java.io.DataOutputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.InputStream;
import java.net.HttpURLConnection;
import java.net.MalformedURLException;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.TreeMap;

import static utilities.Utilities.*;

import solver.GreedyHeuristic;

/**
 *
 * @author yaw
 */
public class DataInOut {

    private String basePath;
    private String dataset;
    private String scenario;
    private DataStorer data;

    public void loadData(String basePath, String dataset, String scenario, DataStorer data) {
        this.basePath = basePath;
        this.dataset = dataset;
        this.scenario = scenario;
        this.data = data;

        System.out.println("Loading Geography...");
        loadGeography();
        System.out.println("Loading Source Data...");
        loadSources();
        System.out.println("Loading Sink Data...");
        loadSinks();
        System.out.println("Loading Transport Data...");
        loadTransport();
        System.out.print("Loading Delaunay Pairs...");
        loadDelaunayPairs();
        System.out.print("Loading Candidate Graph...");
        loadCandidateGraph();
        System.out.println("Data Loaded.");
    }

    private void loadGeography() {
        String path = basePath + "/" + dataset + "/BaseData/CostNetwork/Construction Costs.txt";
        try (BufferedReader br = new BufferedReader(new FileReader(path))) {
            br.readLine();
            br.readLine();

            // Read dimensions.
            String line = br.readLine();
            String[] elements = line.split("\\s+");
            data.setWidth(Integer.parseInt(elements[1]));

            line = br.readLine();
            elements = line.split("\\s+");
            data.setHeight(Integer.parseInt(elements[1]));

            // Read conversions.
            line = br.readLine();
            elements = line.split("\\s+");
            data.setLowerLeftX(Double.parseDouble(elements[1]));

            line = br.readLine();
            elements = line.split("\\s+");
            data.setLowerLeftY(Double.parseDouble(elements[1]));

            line = br.readLine();
            elements = line.split("\\s+");
            data.setCellSize(Double.valueOf(elements[1]));
        } catch (IOException e) {
            path = basePath + "/" + dataset + "/BaseData/CostNetwork/Construction Costs.csv";
            try (BufferedReader br = new BufferedReader(new FileReader(path))) {
                br.readLine();
                br.readLine();

                // Read dimensions.
                String line = br.readLine();
                String[] elements = line.split(",");
                data.setWidth(Integer.parseInt(elements[1]));

                line = br.readLine();
                elements = line.split(",");
                data.setHeight(Integer.parseInt(elements[1]));

                // Read conversions.
                line = br.readLine();
                elements = line.split(",");
                data.setLowerLeftX(Double.parseDouble(elements[1]));

                line = br.readLine();
                elements = line.split(",");
                data.setLowerLeftY(Double.parseDouble(elements[1]));

                line = br.readLine();
                elements = line.split(",");
                data.setCellSize(Double.valueOf(elements[1]));
            } catch (IOException e2) {
                System.out.println(e.getMessage());
            }
        }
    }

    public void loadCosts() {
        double[][] rightOfWayCosts = new double[0][0];
        double[][] constructionCosts = new double[0][0];
        double[][] routingCosts = new double[0][0];

        String path = basePath + "/" + dataset + "/BaseData/CostNetwork/Construction Costs.csv";

        long startTime = System.nanoTime();
        // Load construction costs from csv file.
        try (BufferedReader br = new BufferedReader(new FileReader(path))) {
            System.out.println("Loading " + path);
            // Create construction costs array.
            constructionCosts = new double[data.getWidth() * data.getHeight() + 1][8];
            for (int i = 0; i < constructionCosts.length; i++) {
                for (int j = 0; j < constructionCosts[i].length; j++) {
                    constructionCosts[i][j] = Double.MAX_VALUE;
                }
            }

            for (int i = 0; i < 8; i++) {
                br.readLine();
            }

            String line = br.readLine();
            while (line != null) {
                String costLine = br.readLine();

                int currentCellIndex = 0;
                int nextCellIndex = line.indexOf(",");
                int centerCell = Integer.parseInt(line.substring(currentCellIndex, nextCellIndex));
                currentCellIndex = nextCellIndex + 1;
                int currentCostIndex = 0;
                int nextCostIndex = 0;
                boolean moreNeighbors = true;
                while (moreNeighbors) {
                    nextCellIndex = line.indexOf(",", currentCellIndex);
                    if (nextCellIndex == -1) {
                        moreNeighbors = false;
                        nextCellIndex = line.length();
                    }
                    int neighborCell = Integer.parseInt(line.substring(currentCellIndex, nextCellIndex));
                    currentCellIndex = nextCellIndex + 1;

                    nextCostIndex = costLine.indexOf(",", currentCostIndex);
                    if (nextCostIndex == -1) {
                        nextCostIndex = costLine.length();
                    }
                    double cost = Double.parseDouble(costLine.substring(currentCostIndex, nextCostIndex));
                    currentCostIndex = nextCostIndex + 1;

                    constructionCosts[centerCell][data.getNeighborNum(centerCell, neighborCell)] = cost;
                }

                line = br.readLine();
            }
            long endTime = System.nanoTime();
            System.out.println("Load construction csv file: " + (endTime - startTime)/1e9 + " s");

        } catch (IOException e) {
            // Load construction costs from text file.
            path = basePath + "/" + dataset + "/BaseData/CostNetwork/Construction Costs.txt";
            try (BufferedReader br = new BufferedReader(new FileReader(path))) {
                // Create construction costs array.
                constructionCosts = new double[data.getWidth() * data.getHeight() + 1][8];
                for (int i = 0; i < constructionCosts.length; i++) {
                    for (int j = 0; j < constructionCosts[i].length; j++) {
                        constructionCosts[i][j] = Double.MAX_VALUE;
                    }
                }

                for (int i = 0; i < 8; i++) {
                    br.readLine();
                }

                String line = br.readLine();
                while (line != null) {
                    String costLine = br.readLine();
                    String[] costs = costLine.split("\\s+");
                    String[] cells = line.split("\\s+");

                    int centerCell = Integer.parseInt(cells[0]);
                    for (int i = 1; i < costs.length; i++) {
                        constructionCosts[centerCell][data.getNeighborNum(centerCell, Integer.parseInt(cells[i]))] = Double.parseDouble(costs[i]);
                    }
                    line = br.readLine();
                }
            } catch (IOException e2) {
                System.out.println(e2.getMessage());
            }
        }

        // Load right of way costs.  
        path = basePath + "/" + dataset + "/BaseData/CostNetwork/RightOfWay Costs.txt";
        try (BufferedReader br = new BufferedReader(new FileReader(path))) {
            // Create right of way cost array.
            rightOfWayCosts = new double[data.getWidth() * data.getHeight() + 1][8];
            for (int i = 0; i < rightOfWayCosts.length; i++) {
                for (int j = 0; j < rightOfWayCosts[i].length; j++) {
                    rightOfWayCosts[i][j] = Double.MAX_VALUE;
                }
            }

            for (int i = 0; i < 8; i++) {
                br.readLine();
            }

            String line = br.readLine();
            while (line != null) {
                String costLine = br.readLine();
                String[] costs = costLine.split("\\s+");
                String[] cells = line.split("\\s+");
                int centerCell = Integer.parseInt(cells[0]);
                for (int i = 1; i < costs.length; i++) {
                    rightOfWayCosts[centerCell][data.getNeighborNum(centerCell, Integer.parseInt(cells[i]))] = Double.parseDouble(costs[i]);
                }
                line = br.readLine();
            }
        } catch (IOException e) {
            rightOfWayCosts = null;
        }

        // Load routing costs.
        path = basePath + "/" + dataset + "/BaseData/CostNetwork/Routing Costs.txt";
        startTime = System.nanoTime();
        routingCosts = new double[data.getWidth() * data.getHeight() + 1][8];
        long endTime = System.nanoTime();
        System.out.println("Setup routing costs array: " + (endTime - startTime)/1e9 + " s");
        try (BufferedReader br = new BufferedReader(new FileReader(path))) {
            for (int i = 0; i < routingCosts.length; i++) {
                for (int j = 0; j < routingCosts[i].length; j++) {
                    routingCosts[i][j] = Double.MAX_VALUE;
                }
            }

            for (int i = 0; i < 8; i++) {
                br.readLine();
            }

            String line = br.readLine();
            while (line != null) {
                String costLine = br.readLine();
                String[] costs = costLine.split("\\s+");
                String[] cells = line.split("\\s+");
                int centerCell = Integer.parseInt(cells[0]);
                for (int i = 1; i < costs.length; i++) {
                    routingCosts[centerCell][data.getNeighborNum(centerCell, Integer.parseInt(cells[i]))] = Double.parseDouble(costs[i]);
                }
                line = br.readLine();
            }
        } catch (IOException e) {
            startTime = System.nanoTime();
            for (int i = 0; i < routingCosts.length; i++) {
                for (int j = 0; j < routingCosts[i].length; j++) {
                    routingCosts[i][j] = constructionCosts[i][j];
                    if (rightOfWayCosts != null) {
                        routingCosts[i][j] += rightOfWayCosts[i][j];
                    }
                }
            }
            endTime = System.nanoTime();
            System.out.println("Populate routing costs array: " + (endTime - startTime)/1e9 + " s");
        }

        startTime = System.nanoTime();
        data.setConstructionCosts(constructionCosts);
        data.setRightOfWayCosts(rightOfWayCosts);
        data.setRoutingCosts(routingCosts);
        endTime = System.nanoTime();
        System.out.println("data.setX() method calls: " + (endTime - startTime)/1e9 + " s");
    }

    private void loadSources() {
        String sourcePath = basePath + "/" + dataset + "/Scenarios/" + scenario + "/Sources/Sources.txt";
        try (BufferedReader br = new BufferedReader(new FileReader(sourcePath))) {
            br.readLine();
            String line = br.readLine();
            ArrayList<Source> sources = new ArrayList<>();
            while (line != null) {
                String[] elements = line.split("\\s+");
                Source source = new Source(data);
                source.setLabel(elements[0]);
                source.setCellNum(data.latLonToCell(Double.parseDouble(elements[7]), Double.parseDouble(elements[6])));
                source.setOpeningCost(Double.parseDouble(elements[1]));
                source.setOMCost(Double.parseDouble(elements[2]));
                source.setCaptureCost(Double.parseDouble(elements[3]));
                source.setProductionRates(csvToDoubleArray(elements[4]));
                sources.add(source);
                line = br.readLine();
            }
            data.setSources(sources.toArray(new Source[0]));
        } catch (IOException e) {
            System.out.println(e.getMessage());
        }
    }

    private void loadSinks() {
        String sinkPath = basePath + "/" + dataset + "/Scenarios/" + scenario + "/Sinks/Sinks.txt";
        try (BufferedReader br = new BufferedReader(new FileReader(sinkPath))) {
            br.readLine();
            String line = br.readLine();
            ArrayList<Sink> sinks = new ArrayList<>();
            while (line != null) {
                String[] elements = line.split("\\s+");
                Sink sink = new Sink(data);
                sink.setLabel(elements[0]);
                sink.setCellNum(data.latLonToCell(Double.parseDouble(elements[11]), Double.parseDouble(elements[10])));
                sink.setOpeningCost(Double.parseDouble(elements[3]));
                sink.setOMCost(Double.parseDouble(elements[4]));
                sink.setWellOpeningCost(Double.parseDouble(elements[6]));
                sink.setWellOMCost(Double.parseDouble(elements[7]));
                sink.setInjectionCost(Double.parseDouble(elements[8]));
                sink.setWellCapacity(Double.parseDouble(elements[5]));
                sink.setCapacities(csvToDoubleArray(elements[2]));
                sinks.add(sink);
                line = br.readLine();
            }
            data.setSinks(sinks.toArray(new Sink[0]));
        } catch (IOException e) {
            System.out.println(e.getMessage());
        }
    }

    private void loadTransport() {
        String transportPath = basePath + "/" + dataset + "/Scenarios/" + scenario + "/Transport/Linear.txt";
        try (BufferedReader br = new BufferedReader(new FileReader(transportPath))) {
            br.readLine();
            String line = br.readLine();
            ArrayList<LinearComponent> linearComponents = new ArrayList<>();
            while (line != null) {
                String[] elements = line.split("\\s+");
                LinearComponent linearComponent = new LinearComponent(data);
                linearComponent.setConSlope(Double.parseDouble(elements[1]));
                linearComponent.setConIntercept(Double.parseDouble(elements[2]));
                if (elements.length > 3) {
                    linearComponent.setRowSlope(Double.parseDouble(elements[3]));
                    linearComponent.setRowIntercept(Double.parseDouble(elements[4]));
                }
                linearComponents.add(linearComponent);
                line = br.readLine();
            }

            // Set max pipeline capacities.
            for (int c = 0; c < linearComponents.size(); c++) {
                double maxCap = data.getMaxAnnualCapturable();  // Do not make Double.MaxValue. CPLEX does not do well with infinity here.
                if (c < linearComponents.size() - 1) {
                    double slope1 = linearComponents.get(c).getConSlope() + linearComponents.get(c).getRowSlope();
                    double intercept1 = linearComponents.get(c).getConIntercept() + linearComponents.get(c).getRowIntercept();
                    double slope2 = linearComponents.get(c + 1).getConSlope() + linearComponents.get(c + 1).getRowSlope();
                    double intercept2 = linearComponents.get(c + 1).getConIntercept() + linearComponents.get(c + 1).getRowIntercept();
                    maxCap = (intercept2 - intercept1) / (slope1 - slope2);
                }
                linearComponents.get(c).setMaxCapacity(maxCap);
            }

            data.setLinearComponents(linearComponents.toArray(new LinearComponent[0]));
        } catch (IOException e) {
            System.out.println(e.getMessage());
        }
    }

    private void loadCandidateGraph() {
        // Check if file exists
        String candidateGraphPath = basePath + "/" + dataset + "/Scenarios/" + scenario + "/Network/CandidateNetwork/CandidateNetwork.txt";
        if (new File(candidateGraphPath).exists()) {
            // Load from file.
            try (BufferedReader br = new BufferedReader(new FileReader(candidateGraphPath))) {
                String line = br.readLine();
                // Determine data version
                int routeStarting = 5;
                if (!line.startsWith("Vertex1")) {
                    routeStarting = 4;
                    br.readLine();
                    br.readLine();
                    br.readLine();
                }
                if (!line.contains("ConCost")) {
                    routeStarting = 3;
                }
                line = br.readLine();

                HashSet<Integer> graphVertices = new HashSet<>();
                HashMap<Edge, Double> graphEdgeCosts = new HashMap<>();
                HashMap<Edge, Double> graphEdgeConstructionCosts = new HashMap<>();
                HashMap<Edge, Double> graphEdgeRightOfWayCosts = new HashMap<>();

                HashMap<Edge, int[]> graphEdgeRoutes = new HashMap<>();
                while (line != null) {
                    String[] elements = line.split("\\s+");
                    int v1 = Integer.parseInt(elements[0]);
                    int v2 = Integer.parseInt(elements[1]);
                    Edge edge = new Edge(v1, v2);
                    graphVertices.add(v1);
                    graphVertices.add(v2);
                    double cost = Double.parseDouble(elements[2]);

                    double conCost = 0;
                    double rowCost = 0;
                    if (routeStarting == 5) {
                        conCost = Double.parseDouble(elements[3]);
                        rowCost = Double.parseDouble(elements[4]);
                    }

                    ArrayList<Integer> route = new ArrayList<>();
                    for (int i = routeStarting; i < elements.length; i++) {
                        route.add(Integer.parseInt(elements[i]));
                    }

                    graphEdgeCosts.put(edge, cost);
                    graphEdgeRoutes.put(edge, convertIntegerArray(route.toArray(new Integer[0])));

                    if (routeStarting == 5) {
                        graphEdgeConstructionCosts.put(edge, conCost);
                        graphEdgeRightOfWayCosts.put(edge, rowCost);
                    }

                    // Prepare for next entry
                    line = br.readLine();
                }

                int[] vertices = new int[graphVertices.size()];
                int i = 0;
                for (int vertex : graphVertices) {
                    vertices[i++] = vertex;
                }
                Arrays.sort(vertices);

                data.setGraphVertices(vertices);
                data.setGraphEdgeCosts(graphEdgeCosts);
                data.setGraphEdgeRoutes(graphEdgeRoutes);

                if (routeStarting == 5) {
                    data.setGraphEdgeConstructionCosts(graphEdgeConstructionCosts);
                    data.setGraphEdgeRightOfWayCosts(graphEdgeRightOfWayCosts);
                }

                System.out.println();
            } catch (IOException e) {
                System.out.println(e.getMessage());
            }
        } else {
            System.out.println("Not Yet Generated.");
        }
    }

    private void loadDelaunayPairs() {
        // Check if file exists
        String delaunayPairsPath = basePath + "/" + dataset + "/Scenarios/" + scenario + "/Network/DelaunayNetwork/DelaunayPaths.txt";
        if (new File(delaunayPairsPath).exists()) {
            // Load from file.
            try (BufferedReader br = new BufferedReader(new FileReader(delaunayPairsPath))) {
                br.readLine();
                String line = br.readLine();

                HashSet<Edge> pairs = new HashSet<>();
                while (line != null) {
                    String[] elements = line.split("\\s+");
                    int v1 = Integer.parseInt(elements[4]);
                    int v2 = Integer.parseInt(elements[5]);
                    Edge edge = new Edge(v1, v2);
                    pairs.add(edge);

                    // Prepare for next entry
                    line = br.readLine();
                }

                data.setDelaunayPairs(pairs);
                System.out.println();
            } catch (IOException e) {
                System.out.println(e.getMessage());
            }
        } else {
            System.out.println("Not Yet Generated.");
        }
    }

    public void loadPriceConfiguration() {
        // Check if file exists
        String pricesPath = basePath + "/" + dataset + "/Scenarios/" + scenario + "/Configurations/priceInput.csv";
        if (new File(pricesPath).exists()) {
            // Load from file.
            try (BufferedReader br = new BufferedReader(new FileReader(pricesPath))) {
                br.readLine();
                String line = br.readLine();

                String[] elements = line.split(",");
                double min = Double.parseDouble(elements[0]);
                double max = Double.parseDouble(elements[1]);
                double step = Double.parseDouble(elements[2]);

                // Make prices array
                int num = (int) Math.floor((max - min + 1) / step);
                double[] prices = new double[num];
                for (int i = 0; i < prices.length; i++) {
                    prices[i] = min + i * step;
                }
                data.setPriceConfiguration(prices);
            } catch (IOException e) {
                System.out.println(e.getMessage());
            }
        } else {
            System.out.println("No price configuration file.");
        }
    }

    public void loadTimeConfiguration() {
        // Check if file exists
        String timeConfigurationPath = basePath + "/" + dataset + "/Scenarios/" + scenario + "/Configurations/timeInput.csv";
        if (new File(timeConfigurationPath).exists()) {
            // Load from file.
            try (BufferedReader br = new BufferedReader(new FileReader(timeConfigurationPath))) {
                ArrayList<double[]> timeEntries = new ArrayList<>();
                br.readLine();
                String line = br.readLine();
                while (line != null) {
                    String[] elements = line.split(",");
                    double timeslot = Double.parseDouble(elements[0]);
                    double numYears = Double.parseDouble(elements[1]);
                    double capTarget = Double.parseDouble(elements[2]);

                    timeEntries.add(new double[]{timeslot, numYears, capTarget});
                    line = br.readLine();
                }

                // Store time configuration
                double[][] timeConfiguration = new double[timeEntries.size()][3];
                for (int i = 0; i < timeEntries.size(); i++) {
                    timeConfiguration[i] = timeEntries.get(i);
                }
                data.setTimeConfiguration(timeConfiguration);
            } catch (IOException e) {
                System.out.println(e.getMessage());
            }
        } else {
            System.out.println("No time configuration file.");
        }
    }

    public void saveDelaunayPairs() {
        HashSet<Edge> delaunayPairs = data.getDelaunayPairs();

        String delaunayPairsPath = basePath + "/" + dataset + "/Scenarios/" + scenario + "/Network/DelaunayNetwork/DelaunayPaths.txt";

        // Save to file.
        try (BufferedWriter bw = new BufferedWriter(new FileWriter(delaunayPairsPath))) {
            bw.write("#  Selected node pairs\n");
            for (Edge pair : delaunayPairs) {
                int vNum = data.sourceNum(pair.v1);
                if (vNum > -1) {
                    bw.write("SOURCE\t" + data.getSources()[vNum].getLabel() + "\t");
                } else {
                    bw.write("SINK\t" + data.getSinks()[data.sinkNum(pair.v1)].getLabel() + "\t");
                }
                vNum = data.sourceNum(pair.v2);
                if (vNum > -1) {
                    bw.write("SOURCE\t" + data.getSources()[vNum].getLabel() + "\t");
                } else {
                    bw.write("SINK\t" + data.getSinks()[data.sinkNum(pair.v2)].getLabel() + "\t");
                }
                bw.write(pair.v1 + "\t" + pair.v2 + "\n");
            }
        } catch (IOException e) {
            System.out.println(e.getMessage());
        }
    }

    public void saveCandidateGraph() {
        HashMap<Edge, Double> graphEdgeCosts = data.getGraphEdgeCosts();
        HashMap<Edge, int[]> graphEdgeRoutes = data.getGraphEdgeRoutes();
        HashMap<Edge, Double> graphEdgeConstructionCosts = data.getGraphEdgeConstructionCosts();
        HashMap<Edge, Double> graphEdgeRightOfWayCosts = data.getGraphEdgeRightOfWayCosts();

        String rawPathsPath = basePath + "/" + dataset + "/Scenarios/" + scenario + "/Network/CandidateNetwork/CandidateNetwork.txt";

        // Save to file.
        try (BufferedWriter bw = new BufferedWriter(new FileWriter(rawPathsPath))) {
            bw.write("Vertex1\tVertex2\tCost\tConCost\tROWCost\tCellRoute\n");
            for (Edge e : graphEdgeRoutes.keySet()) {
                bw.write(e.v1 + "\t" + e.v2 + "\t" + graphEdgeCosts.get(e) + "\t" + graphEdgeConstructionCosts.get(e) + "\t" + graphEdgeRightOfWayCosts.get(e));
                int[] route = graphEdgeRoutes.get(e);
                for (int vertex : route) {
                    bw.write("\t" + vertex);
                }
                bw.write("\n");
            }
        } catch (IOException e) {
            System.out.println(e.getMessage());
        }
    }

    // Heuristic
    public void saveHeuristicSolution(File solutionDirectory, GreedyHeuristic heuristic) {
        // Collect data.
        Source[] sources = heuristic.getSources();
        Sink[] sinks = heuristic.getSinks();
        int[] graphVertices = heuristic.getGraphVertices();
        HeuristicEdge[][] adjacencyMatrix = heuristic.getAdjacencyMatrix();
        HashMap<Integer, Integer> cellNumToVertexNum = heuristic.getCellVertexMap();
        double crf = data.getCrf();

        try (BufferedWriter bw = new BufferedWriter(new FileWriter(solutionDirectory.toString() + "/solution.txt"))) {
            bw.write("crf:\t" + crf + "\n");
            bw.write("captureTarget:\t" + data.getTargetCaptureAmount() + "\n");
            bw.write("projectLength:\t" + data.getProjectLength() + "\n");

            bw.write("SourceCell\tSourceLabel\tCaptureAmount\tCost\n");
            for (Source src : sources) {
                if (src.getRemainingCapacity() < src.getProductionRate()) {
                    double captureAmount = src.getProductionRate() - src.getRemainingCapacity();
                    double cost = src.getOpeningCost(crf) + src.getCaptureCost() * captureAmount;
                    bw.write(src.getCellNum() + "\t" + src.getLabel() + "\t" + captureAmount + "\t" + cost + "\n");
                }
            }

            bw.write("Sink\tSinkLabel\tInjectAmount\tCost\n");
            for (Sink snk : sinks) {
                if (snk.getRemainingCapacity() < (snk.getCapacity() / data.getProjectLength())) {
                    double injectAmount = (snk.getCapacity() / data.getProjectLength()) - snk.getRemainingCapacity();
                    double cost = snk.getOpeningCost(crf) + snk.getInjectionCost() * injectAmount + snk.getNumWells() * snk.getWellOpeningCost(crf);
                    bw.write(snk.getCellNum() + "\t" + snk.getLabel() + "\t" + injectAmount + "\t" + cost + "\n");
                }
            }

            bw.write("EdgeSrc\tEdgeSnk\tFlowAmount\tCost\n");
            for (int u = 0; u < graphVertices.length; u++) {
                for (int v = 0; v < graphVertices.length; v++) {
                    HeuristicEdge edge = adjacencyMatrix[u][v];
                    if (edge != null && edge.currentHostingAmount > 0) {
                        double flowAmount = edge.currentHostingAmount;
                        double cost = edge.buildCost[edge.currentSize] + edge.transportCost[edge.currentSize] * flowAmount;
                        bw.write(edge.v1 + "\t" + edge.v2 + "\t" + flowAmount + "\t" + cost + "\n");
                    }
                }
            }
        } catch (IOException e) {
            System.out.println(e.getMessage());
        }
    }

    // Heuristic
    public Solution loadGreedyHeuristicSolution(String solutionPath) {
        Solution soln = new Solution();
        Source[] sources = data.getSources();
        Sink[] sinks = data.getSinks();

        try (BufferedReader br = new BufferedReader(new FileReader(solutionPath.toString() + "/solution.txt"))) {
            String line = br.readLine();
            soln.setCRF(Double.parseDouble(line.split("\t")[1]));

            line = br.readLine();

            line = br.readLine();
            soln.setProjectLength(Integer.parseInt(line.split("\t")[1]));

            line = br.readLine();
            line = br.readLine();
            while (!line.startsWith("Sink")) {
                String[] sourceComponents = line.split("\t");
                Source source = sources[data.sourceNum(Integer.parseInt(sourceComponents[0]))];
                double captureAmount = Double.parseDouble(sourceComponents[2]);
                double cost = Double.parseDouble(sourceComponents[3]);
                soln.addSourceCaptureAmount(source, captureAmount);
                soln.addSourceCostComponent(source, cost);
                line = br.readLine();
            }

            line = br.readLine();
            while (!line.startsWith("EdgeSrc")) {
                String[] sinkComponents = line.split("\t");
                Sink sink = sinks[data.sinkNum(Integer.parseInt(sinkComponents[0]))];
                double injectAmount = Double.parseDouble(sinkComponents[2]);
                double cost = Double.parseDouble(sinkComponents[3]);
                soln.addSinkStorageAmount(sink, injectAmount);
                soln.addSinkCostComponent(sink, cost);
                line = br.readLine();
            }

            line = br.readLine();
            while (line != null) {
                String[] edgeComponents = line.split("\t");
                Edge edge = new Edge(Integer.parseInt(edgeComponents[0]), Integer.parseInt(edgeComponents[1]));
                double flowAmount = Double.parseDouble(edgeComponents[2]);
                double cost = Double.parseDouble(edgeComponents[3]);
                soln.addEdgeTransportAmount(edge, flowAmount);
                soln.addEdgeCostComponent(edge, cost);
                line = br.readLine();
            }
        } catch (IOException e) {
            System.out.println(e.getMessage());
        }

        return soln;
    }

    public Solution loadSolution(String solutionPath, int timeslot) {
        double threshold = .000001;
        Solution soln = new Solution();

        boolean ilp = true; //ILP vs LP

        // Make file paths.
        File solFile = null;
        File mpsFile = null;
        for (File f : new File(solutionPath).listFiles()) {
            if (f.getName().endsWith(".sol")) {
                solFile = f;
            } else if (f.getName().endsWith(".mps")) {
                mpsFile = f;
                if (f.getName().startsWith("flowCap")) {
                    ilp = false;
                }
            }
        }

        // Collect data.
        Source[] sources = data.getSources();
        Sink[] sinks = data.getSinks();
        int[] graphVertices = data.getGraphVertices();
        HashMap<Edge, Double> edgeConstructionCosts = data.getGraphEdgeConstructionCosts();

        // Make cell/index maps.
        HashMap<Source, Integer> sourceCellToIndex = new HashMap<>();
        HashMap<Integer, Source> sourceIndexToCell = new HashMap<>();
        HashMap<Sink, Integer> sinkCellToIndex = new HashMap<>();
        HashMap<Integer, Sink> sinkIndexToCell = new HashMap<>();
        HashMap<Integer, Integer> vertexCellToIndex = new HashMap<>();
        HashMap<Integer, Integer> vertexIndexToCell = new HashMap<>();
        HashMap<UnidirEdge, Integer> edgeToIndex = new HashMap<>();
        HashMap<Integer, UnidirEdge> edgeIndexToEdge = new HashMap<>();

        // Initialize cell/index maps.
        for (int i = 0; i < sources.length; i++) {
            sourceCellToIndex.put(sources[i], i);
            sourceIndexToCell.put(i, sources[i]);
        }
        for (int i = 0; i < sinks.length; i++) {
            sinkCellToIndex.put(sinks[i], i);
            sinkIndexToCell.put(i, sinks[i]);
        }
        for (int i = 0; i < graphVertices.length; i++) {
            vertexCellToIndex.put(graphVertices[i], i);
            vertexIndexToCell.put(i, graphVertices[i]);
        }
        int index = 0;
        for (Edge e : edgeConstructionCosts.keySet()) {
            UnidirEdge e1 = new UnidirEdge(e.v1, e.v2);
            edgeToIndex.put(e1, index);
            edgeIndexToEdge.put(index, e1);
            index++;

            UnidirEdge e2 = new UnidirEdge(e.v2, e.v1);
            edgeToIndex.put(e2, index);
            edgeIndexToEdge.put(index, e2);
            index++;
        }

        HashMap<String, Double> variableValues = new HashMap<>();
        HashMap<Integer, Integer> timeslotLengths = new HashMap<>();

        try (BufferedReader br = new BufferedReader(new FileReader(solFile))) {
            String line = br.readLine();
            while (!line.equals(" <variables>")) {
                line = br.readLine();
            }
            line = br.readLine();

            while (!line.equals(" </variables>")) {
                String[] partition = line.split("\"");
                String[] variable;
                if (ilp) {
                    variable = new String[]{partition[1], partition[3], partition[5]};
                } else {
                    variable = new String[]{partition[1], partition[3], partition[7]};
                }

                if (Double.parseDouble(variable[2]) > threshold) {
                    variableValues.put(variable[0], Double.parseDouble(variable[2]));
                    String[] components = variable[0].split("\\]\\[|\\[|\\]");
                    if (timeslot == -1) {
                        if (components[0].equals("a")) {
                            soln.addSourceCaptureAmount(sources[Integer.parseInt(components[1])], Double.parseDouble(variable[2]));
                        } else if (components[0].equals("b")) {
                            soln.addSinkStorageAmount(sinks[Integer.parseInt(components[1])], Double.parseDouble(variable[2]));
                        } else if (components[0].equals("p")) {
                            if (components.length == 4) {
                                soln.addEdgeTransportAmount(new Edge(vertexIndexToCell.get(Integer.parseInt(components[1])), vertexIndexToCell.get(Integer.parseInt(components[2]))), Double.parseDouble(variable[2]));  
                                soln.setEdgeTrend(new Edge(vertexIndexToCell.get(Integer.parseInt(components[1])), vertexIndexToCell.get(Integer.parseInt(components[2]))), Integer.parseInt(components[3]));
                            } else {
                                UnidirEdge unidirEdge = edgeIndexToEdge.get(Integer.parseInt(components[1]));
                                soln.addEdgeTransportAmount(new Edge(unidirEdge.v1, unidirEdge.v2), Double.parseDouble(variable[2]));
                                soln.setEdgeTrend(new Edge(unidirEdge.v1, unidirEdge.v2), Integer.parseInt(components[2]));
                            }
                        } else if (components[0].equals("w")) {
                            soln.addSinkNumWells(sinks[Integer.parseInt(components[1])], (int) Math.round(Double.parseDouble(variable[2])));
                        }
                    } else {
                        if (components[0].equals("a") && (Integer.parseInt(components[2]) == timeslot)) {
                            soln.addSourceCaptureAmount(sources[Integer.parseInt(components[1])], Double.parseDouble(variable[2]));
                        } else if (components[0].equals("b") && (Integer.parseInt(components[2]) == timeslot)) {
                            soln.addSinkStorageAmount(sinks[Integer.parseInt(components[1])], Double.parseDouble(variable[2]));
                        } else if (components[0].equals("x")) {
                            if (components.length == 5) {
                                if (Integer.parseInt(components[4]) == timeslot) {
                                    soln.addEdgeTransportAmount(new Edge(vertexIndexToCell.get(Integer.parseInt(components[1])), vertexIndexToCell.get(Integer.parseInt(components[2]))), Double.parseDouble(variable[2]));
                                }
                            } else {
                                if (Integer.parseInt(components[3]) == timeslot) {
                                    UnidirEdge unidirEdge = edgeIndexToEdge.get(Integer.parseInt(components[1]));
                                    soln.addEdgeTransportAmount(new Edge(unidirEdge.v1, unidirEdge.v2), Double.parseDouble(variable[2]));
                                }
                            }
                        } else if (components[0].equals("w") && (Integer.parseInt(components[2]) == timeslot)) {
                            soln.addSinkNumWells(sinks[Integer.parseInt(components[1])], Integer.parseInt(variable[2]));
                        }
                    }

                    if (variable[0].equals("crf")) {
                        soln.setCRF(Double.parseDouble(variable[2]));
                    } else if (variable[0].equals("projectLength")) {
                        soln.setProjectLength(Integer.parseInt(variable[2]));
                    } else if (variable[0].startsWith("N")) {
                        int t = Integer.parseInt(variable[0].substring(1));
                        int length = Integer.parseInt(variable[2]);
                        timeslotLengths.put(t, length);
                        if (t == timeslot) {
                            soln.setProjectLength(length);
                        }
                    }
                }
                line = br.readLine();
            }
        } catch (IOException e) {
            System.out.println(e.getMessage());
        }

        // load costs into solution.
        if (ilp) {
            try (BufferedReader br = new BufferedReader(new FileReader(mpsFile))) {
                String line = br.readLine();
                while (!line.equals("COLUMNS")) {
                    line = br.readLine();
                }
                br.readLine();
                line = br.readLine();

                while (!line.equals("RHS")) {
                    String[] column = line.replaceFirst("\\s+", "").split("\\s+");
                    if (column[1].equals("OBJ") && variableValues.keySet().contains(column[0])) {
                        String[] components = column[0].split("\\]\\[|\\[|\\]");
                        if (timeslot == -1) {
                            if (column[0].charAt(0) == 's' || column[0].charAt(0) == 'a') {
                                double cost = variableValues.get(column[0]) * Double.parseDouble(column[2]);
                                soln.addSourceCostComponent(sources[Integer.parseInt(components[1])], cost);
                            } else if (column[0].charAt(0) == 'r' || column[0].charAt(0) == 'w' || column[0].charAt(0) == 'b') {
                                double cost = variableValues.get(column[0]) * Double.parseDouble(column[2]);
                                soln.addSinkCostComponent(sinks[Integer.parseInt(components[1])], cost);
                            } else if (column[0].charAt(0) == 'p' || column[0].charAt(0) == 'y') {
                                double cost = variableValues.get(column[0]) * Double.parseDouble(column[2]);
                                if (components.length == 4) {
                                    soln.addEdgeCostComponent(new Edge(vertexIndexToCell.get(Integer.parseInt(components[1])), vertexIndexToCell.get(Integer.parseInt(components[2]))), cost);
                                } else {
                                    UnidirEdge unidirEdge = edgeIndexToEdge.get(Integer.parseInt(components[1]));
                                    soln.addEdgeCostComponent(new Edge(unidirEdge.v1, unidirEdge.v2), cost);
                                }
                            }
                        } else {
                            if ((column[0].charAt(0) == 's' || column[0].charAt(0) == 'a') && (Integer.parseInt(components[2]) == timeslot)) {
                                double cost = variableValues.get(column[0]) * Double.parseDouble(column[2]) / timeslotLengths.get(timeslot);
                                soln.addSourceCostComponent(sources[Integer.parseInt(components[1])], cost);
                            } else if ((column[0].charAt(0) == 'r' || column[0].charAt(0) == 'w' || column[0].charAt(0) == 'b') && (Integer.parseInt(components[2]) == timeslot)) {
                                double cost = variableValues.get(column[0]) * Double.parseDouble(column[2]) / timeslotLengths.get(timeslot);
                                soln.addSinkCostComponent(sinks[Integer.parseInt(components[1])], cost);
                            } else if (column[0].charAt(0) == 'p' || column[0].charAt(0) == 'y') {
                                if (components.length == 5) {
                                    // Need to account for pipelines still being paid off in current timeslot.
                                    int timeslotOpened = Integer.parseInt(components[4]);
                                    if (timeslotOpened <= timeslot) {
                                        // Find time remaing.
                                        int timeRemaining = 0;
                                        for (int t = timeslotOpened; t < timeslotLengths.size(); t++) {
                                            timeRemaining += timeslotLengths.get(t);
                                        }
                                        double cost = variableValues.get(column[0]) * Double.parseDouble(column[2]) / timeRemaining;
                                        soln.addEdgeCostComponent(new Edge(vertexIndexToCell.get(Integer.parseInt(components[1])), vertexIndexToCell.get(Integer.parseInt(components[2]))), cost);
                                    }
                                } else {
                                    int timeslotOpened = Integer.parseInt(components[3]);
                                    if (timeslotOpened <= timeslot) {
                                        int timeRemaining = 0;
                                        for (int t = timeslotOpened; t < timeslotLengths.size(); t++) {
                                            timeRemaining += timeslotLengths.get(t);
                                        }
                                        UnidirEdge unidirEdge = edgeIndexToEdge.get(Integer.parseInt(components[1]));
                                        double cost = variableValues.get(column[0]) * Double.parseDouble(column[2]) / timeRemaining;
                                        soln.addEdgeCostComponent(new Edge(unidirEdge.v1, unidirEdge.v2), cost);
                                    }
                                }
                            }
                        }
                    }
                    line = br.readLine();
                }
            } catch (IOException e) {
                System.out.println(e.getMessage());
            }
        } else {
            soln.setSolutionCosts(data);
        }

        return soln;
    }

    public int determineNumTimeslots(String mpsFilePath) {
        File mpsFile = new File(mpsFilePath);
        HashSet<Integer> timeslots = new HashSet<>();

        try (BufferedReader br = new BufferedReader(new FileReader(mpsFile))) {
            String line = br.readLine();
            while (!line.equals("COLUMNS")) {
                line = br.readLine();
            }
            br.readLine();
            line = br.readLine();

            while (!line.equals("RHS")) {
                String[] column = line.replaceFirst("\\s+", "").split("\\s+");
                if (column[0].charAt(0) == 'a') {
                    String[] components = column[0].split("\\]\\[|\\[|\\]");
                    timeslots.add(Integer.parseInt(components[2]));
                }
                line = br.readLine();
            }
        } catch (IOException e) {
            System.out.println(e.getMessage());
        }

        return timeslots.size();
    }

    public void makeShapeFiles(String path, Solution soln) {
        // Make shapefiles if they do not already exist.
        File newDir = new File(path + "/shapeFiles/");
        if (!newDir.exists()) {
            newDir.mkdir();

            // Collect data.
            Source[] sources = data.getSources();
            Sink[] sinks = data.getSinks();
            HashMap<Source, Double> sourceCaptureAmounts = soln.getSourceCaptureAmounts();
            HashMap<Sink, Double> sinkStorageAmounts = soln.getSinkStorageAmounts();
            HashMap<Edge, Double> edgeTransportAmounts = soln.getEdgeTransportAmounts();
            HashMap<Edge, int[]> graphEdgeRoutes = data.getGraphEdgeRoutes();

            // Make source shapefiles.
            EsriPointList sourceList = new EsriPointList();
            String[] sourceAttributeNames = {"Id", "X", "Y", "CO2Cptrd", "MxSpply", "PieWdge", "GensUsed", "MaxGens", "ActlCst", "TtlCst", "Name", "Cell#"};
            int[] sourceAttributeDecimals = {0, 6, 6, 6, 6, 6, 0, 0, 0, 0, 0, 0};
            DbfTableModel sourceAttributeTable = new DbfTableModel(sourceAttributeNames.length);   //12
            for (int colNum = 0; colNum < sourceAttributeNames.length; colNum++) {
                sourceAttributeTable.setColumnName(colNum, sourceAttributeNames[colNum]);
                sourceAttributeTable.setDecimalCount(colNum, (byte) sourceAttributeDecimals[colNum]);
                sourceAttributeTable.setLength(colNum, 10);
                if (sourceAttributeNames[colNum].equals("Id")) {
                    sourceAttributeTable.setType(colNum, DbfTableModel.TYPE_CHARACTER);
                } else {
                    sourceAttributeTable.setType(colNum, DbfTableModel.TYPE_NUMERIC);
                }
            }

            // Order sources.
            TreeMap<Double, ArrayList<Source>> orderedSources = new TreeMap<>();
            for (Source src : sources) {
                if (orderedSources.get(-src.getProductionRate()) == null) {
                    orderedSources.put(-src.getProductionRate(), new ArrayList<Source>());
                }
                orderedSources.get(-src.getProductionRate()).add(src);
            }

            for (ArrayList<Source> sameCapSources : orderedSources.values()) {
                for (Source src : sameCapSources) {
                    EsriPoint source = new EsriPoint(data.cellToLatLon(src.getCellNum())[0], data.cellToLatLon(src.getCellNum())[1]);
                    sourceList.add(source);

                    // Add attributes.
                    ArrayList row = new ArrayList();
                    row.add(src.getLabel());
                    row.add(data.cellToLatLon(src.getCellNum())[1]);
                    row.add(data.cellToLatLon(src.getCellNum())[0]);
                    if (sourceCaptureAmounts.containsKey(src)) {
                        row.add(sourceCaptureAmounts.get(src));
                        row.add(src.getProductionRate());
                        row.add(src.getProductionRate() - sourceCaptureAmounts.get(src));
                    } else {
                        row.add(0);
                        row.add(src.getProductionRate());
                        row.add(src.getProductionRate());
                    }
                    for (int i = 0; i < 6; i++) {
                        row.add(0);
                    }

                    sourceAttributeTable.addRecord(row);
                }
            }

            EsriShapeExport writeSourceShapefiles = new EsriShapeExport(sourceList, sourceAttributeTable, newDir.toString() + "/Sources");
            writeSourceShapefiles.export();
            makeProjectionFile("Sources", newDir.toString());

            // Make sink shapefiles.
            EsriPointList sinkList = new EsriPointList();
            String[] sinkAttributeNames = {"Id", "X", "Y", "CO2Strd", "MxStrg", "PieWdge", "WllsUsed", "MxWlls", "ActCst", "TtlCst", "Name", "Cell#"};
            int[] sinkAttributeDecimals = {0, 6, 6, 6, 6, 6, 0, 0, 0, 0, 0, 0};
            DbfTableModel sinkAttributeTable = new DbfTableModel(sinkAttributeNames.length);   //12
            for (int colNum = 0; colNum < sinkAttributeNames.length; colNum++) {
                sinkAttributeTable.setColumnName(colNum, sinkAttributeNames[colNum]);
                sinkAttributeTable.setDecimalCount(colNum, (byte) sinkAttributeDecimals[colNum]);
                sinkAttributeTable.setLength(colNum, 10);
                if (sinkAttributeNames[colNum].equals("Id")) {
                    sinkAttributeTable.setType(colNum, DbfTableModel.TYPE_CHARACTER);
                } else {
                    sinkAttributeTable.setType(colNum, DbfTableModel.TYPE_NUMERIC);
                }
            }

            // Order sinks.
            TreeMap<Double, ArrayList<Sink>> orderedSinks = new TreeMap<>();
            for (Sink snk : sinks) {
                if (orderedSinks.get(-snk.getCapacity()) == null) {
                    orderedSinks.put(-snk.getCapacity(), new ArrayList<Sink>());
                }
                orderedSinks.get(-snk.getCapacity()).add(snk);
            }

            for (ArrayList<Sink> sameCapSinks : orderedSinks.values()) {
                for (Sink snk : sameCapSinks) {
                    EsriPoint source = new EsriPoint(data.cellToLatLon(snk.getCellNum())[0], data.cellToLatLon(snk.getCellNum())[1]);
                    sinkList.add(source);

                    // Add attributes.
                    ArrayList row = new ArrayList();
                    row.add(snk.getLabel());
                    row.add(data.cellToLatLon(snk.getCellNum())[1]);
                    row.add(data.cellToLatLon(snk.getCellNum())[0]);
                    if (sinkStorageAmounts.containsKey(snk)) {
                        row.add(sinkStorageAmounts.get(snk));
                        row.add(snk.getCapacity() / soln.getProjectLength());
                        row.add(snk.getCapacity() / soln.getProjectLength() - sinkStorageAmounts.get(snk));
                    } else {
                        row.add(0);
                        row.add(snk.getCapacity() / soln.getProjectLength());
                        row.add(snk.getCapacity() / soln.getProjectLength());
                    }
                    for (int i = 0; i < 6; i++) {
                        row.add(0);
                    }

                    sinkAttributeTable.addRecord(row);
                }
            }

            EsriShapeExport writeSinkShapefiles = new EsriShapeExport(sinkList, sinkAttributeTable, newDir.toString() + "/Sinks");
            writeSinkShapefiles.export();
            makeProjectionFile("Sinks", newDir.toString());

            // Make network shapefiles.
            EsriPolylineList edgeList = new EsriPolylineList();
            String[] edgeAttributeNames = {"Id", "CapID", "CapValue", "Flow", "Cost", "LengKM", "LengROW", "LengCONS", "Variable"};
            int[] edgeAttributeDecimals = {0, 0, 0, 6, 0, 0, 0, 0, 0};
            DbfTableModel edgeAttributeTable = new DbfTableModel(edgeAttributeNames.length);   //12
            for (int colNum = 0; colNum < edgeAttributeNames.length; colNum++) {
                edgeAttributeTable.setColumnName(colNum, edgeAttributeNames[colNum]);
                edgeAttributeTable.setDecimalCount(colNum, (byte) edgeAttributeDecimals[colNum]);
                edgeAttributeTable.setLength(colNum, 10);
                if (edgeAttributeNames[colNum].equals("Id")) {
                    edgeAttributeTable.setType(colNum, DbfTableModel.TYPE_CHARACTER);
                } else {
                    edgeAttributeTable.setType(colNum, DbfTableModel.TYPE_NUMERIC);
                }
            }
            for (Edge edg : soln.getOpenedEdges()) {
                // Build route
                int[] route = graphEdgeRoutes.get(edg);
                double[] routeLatLon = new double[route.length * 2];    // Route cells translated into: lat, lon, lat, lon,...
                for (int i = 0; i < route.length; i++) {
                    int cell = route[i];
                    routeLatLon[i * 2] = data.cellToLatLon(cell)[0];
                    routeLatLon[i * 2 + 1] = data.cellToLatLon(cell)[1];
                }

                EsriPolyline edge = new EsriPolyline(routeLatLon, OMGraphic.DECIMAL_DEGREES, OMGraphic.LINETYPE_STRAIGHT);
                edgeList.add(edge);

                // Add attributes.
                ArrayList row = new ArrayList();
                for (int i = 0; i < 3; i++) {
                    row.add(0);
                }
                row.add(edgeTransportAmounts.get(edg));
                for (int i = 0; i < 5; i++) {
                    row.add(0);
                }

                edgeAttributeTable.addRecord(row);
            }

            EsriShapeExport writeEdgeShapefiles = new EsriShapeExport(edgeList, edgeAttributeTable, newDir.toString() + "/Network");
            writeEdgeShapefiles.export();
            makeProjectionFile("Network", newDir.toString());
        }
    }

    public void makeCandidateShapeFiles(String path) {
        // Make shapefiles if they do not already exist.
        File newDir = new File(path + "/shapeFiles/");
        if (!newDir.exists()) {
            newDir.mkdir();

            // Collect data.
            Source[] sources = data.getSources();
            Sink[] sinks = data.getSinks();
            HashMap<Edge, int[]> graphEdgeRoutes = data.getGraphEdgeRoutes();

            // Make source shapefiles.
            EsriPointList sourceList = new EsriPointList();
            String[] sourceAttributeNames = {"Id", "X", "Y"};
            int[] sourceAttributeDecimals = {0, 6, 6};
            DbfTableModel sourceAttributeTable = new DbfTableModel(sourceAttributeNames.length);   //12
            for (int colNum = 0; colNum < sourceAttributeNames.length; colNum++) {
                sourceAttributeTable.setColumnName(colNum, sourceAttributeNames[colNum]);
                sourceAttributeTable.setDecimalCount(colNum, (byte) sourceAttributeDecimals[colNum]);
                sourceAttributeTable.setLength(colNum, 10);
                if (sourceAttributeNames[colNum].equals("Id")) {
                    sourceAttributeTable.setType(colNum, DbfTableModel.TYPE_CHARACTER);
                } else {
                    sourceAttributeTable.setType(colNum, DbfTableModel.TYPE_NUMERIC);
                }
            }

            // Order sources.
            TreeMap<Double, ArrayList<Source>> orderedSources = new TreeMap<>();
            for (Source src : sources) {
                if (orderedSources.get(-src.getProductionRate()) == null) {
                    orderedSources.put(-src.getProductionRate(), new ArrayList<Source>());
                }
                orderedSources.get(-src.getProductionRate()).add(src);
            }

            for (ArrayList<Source> sameCapSources : orderedSources.values()) {
                for (Source src : sameCapSources) {
                    EsriPoint source = new EsriPoint(data.cellToLatLon(src.getCellNum())[0], data.cellToLatLon(src.getCellNum())[1]);
                    sourceList.add(source);

                    // Add attributes.
                    ArrayList row = new ArrayList();
                    row.add(src.getLabel());
                    row.add(data.cellToLatLon(src.getCellNum())[1]);
                    row.add(data.cellToLatLon(src.getCellNum())[0]);

                    sourceAttributeTable.addRecord(row);
                }
            }

            EsriShapeExport writeSourceShapefiles = new EsriShapeExport(sourceList, sourceAttributeTable, newDir.toString() + "/Sources");
            writeSourceShapefiles.export();
            makeProjectionFile("Sources", newDir.toString());

            // Make sink shapefiles.
            EsriPointList sinkList = new EsriPointList();
            String[] sinkAttributeNames = {"Id", "X", "Y"};
            int[] sinkAttributeDecimals = {0, 6, 6};
            DbfTableModel sinkAttributeTable = new DbfTableModel(sinkAttributeNames.length);   //12
            for (int colNum = 0; colNum < sinkAttributeNames.length; colNum++) {
                sinkAttributeTable.setColumnName(colNum, sinkAttributeNames[colNum]);
                sinkAttributeTable.setDecimalCount(colNum, (byte) sinkAttributeDecimals[colNum]);
                sinkAttributeTable.setLength(colNum, 10);
                if (sinkAttributeNames[colNum].equals("Id")) {
                    sinkAttributeTable.setType(colNum, DbfTableModel.TYPE_CHARACTER);
                } else {
                    sinkAttributeTable.setType(colNum, DbfTableModel.TYPE_NUMERIC);
                }
            }

            // Order sinks.
            TreeMap<Double, ArrayList<Sink>> orderedSinks = new TreeMap<>();
            for (Sink snk : sinks) {
                if (orderedSinks.get(-snk.getCapacity()) == null) {
                    orderedSinks.put(-snk.getCapacity(), new ArrayList<Sink>());
                }
                orderedSinks.get(-snk.getCapacity()).add(snk);
            }

            for (ArrayList<Sink> sameCapSinks : orderedSinks.values()) {
                for (Sink snk : sameCapSinks) {
                    EsriPoint source = new EsriPoint(data.cellToLatLon(snk.getCellNum())[0], data.cellToLatLon(snk.getCellNum())[1]);
                    sinkList.add(source);

                    // Add attributes.
                    ArrayList row = new ArrayList();
                    row.add(snk.getLabel());
                    row.add(data.cellToLatLon(snk.getCellNum())[1]);
                    row.add(data.cellToLatLon(snk.getCellNum())[0]);

                    sinkAttributeTable.addRecord(row);
                }
            }

            EsriShapeExport writeSinkShapefiles = new EsriShapeExport(sinkList, sinkAttributeTable, newDir.toString() + "/Sinks");
            writeSinkShapefiles.export();
            makeProjectionFile("Sinks", newDir.toString());

            // Make network shapefiles.
            EsriPolylineList edgeList = new EsriPolylineList();
            String[] edgeAttributeNames = {"Id"};
            int[] edgeAttributeDecimals = {0};
            DbfTableModel edgeAttributeTable = new DbfTableModel(edgeAttributeNames.length);   //12
            for (int colNum = 0; colNum < edgeAttributeNames.length; colNum++) {
                edgeAttributeTable.setColumnName(colNum, edgeAttributeNames[colNum]);
                edgeAttributeTable.setDecimalCount(colNum, (byte) edgeAttributeDecimals[colNum]);
                edgeAttributeTable.setLength(colNum, 10);
                if (edgeAttributeNames[colNum].equals("Id")) {
                    edgeAttributeTable.setType(colNum, DbfTableModel.TYPE_CHARACTER);
                } else {
                    edgeAttributeTable.setType(colNum, DbfTableModel.TYPE_NUMERIC);
                }
            }
            for (Edge edg : graphEdgeRoutes.keySet()) {
                // Build route
                int[] route = graphEdgeRoutes.get(edg);
                double[] routeLatLon = new double[route.length * 2];    // Route cells translated into: lat, lon, lat, lon,...
                for (int i = 0; i < route.length; i++) {
                    int cell = route[i];
                    routeLatLon[i * 2] = data.cellToLatLon(cell)[0];
                    routeLatLon[i * 2 + 1] = data.cellToLatLon(cell)[1];
                }

                EsriPolyline edge = new EsriPolyline(routeLatLon, OMGraphic.DECIMAL_DEGREES, OMGraphic.LINETYPE_STRAIGHT);
                edgeList.add(edge);

                // Add attributes.
                ArrayList row = new ArrayList();
                for (int i = 0; i < 1; i++) {
                    row.add(0);
                }
                edgeAttributeTable.addRecord(row);
            }

            EsriShapeExport writeEdgeShapefiles = new EsriShapeExport(edgeList, edgeAttributeTable, newDir.toString() + "/Network");
            writeEdgeShapefiles.export();
            makeProjectionFile("Network", newDir.toString());
        }
    }

    public void makeProjectionFile(String name, String path) {
        try (BufferedWriter bw = new BufferedWriter(new FileWriter(new File(path, name + ".prj")))) {
            bw.write("GEOGCS[\"GCS_North_American_1983\",DATUM[\"D_North_American_1983\",SPHEROID[\"GRS_1980\",6378137.0,298.257222101]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]]");
        } catch (IOException e) {
            System.out.println(e.getMessage());
        }

    }

    public void makeSolutionFile(String path, Solution soln) {
        try (BufferedWriter bw = new BufferedWriter(new FileWriter(new File(path, "solution.csv")))) {
            bw.write("Project Length," + soln.getProjectLength() + "\n");
            bw.write("CRF," + soln.getCRF() + "\n");
            bw.write("Annual Capture Amount (MTCO2/yr)," + soln.getAnnualCaptureAmount() + "\n");
            bw.write("Total Cost ($M/yr)," + soln.getTotalCost() + "\n");
            bw.write("Capture Cost ($M/yr)," + soln.getTotalAnnualCaptureCost() + "\n");
            bw.write("Transport Cost ($M/yr)," + soln.getTotalAnnualTransportCost() + "\n");
            bw.write("Storage Cost ($M/yr)," + soln.getTotalAnnualStorageCost() + "\n\n");
            bw.write("Source,Capture Amount (MTCO2/yr),Capture Cost ($M/yr)\n");
            HashMap<Source, Double> sourceCaptureAmounts = soln.getSourceCaptureAmounts();
            HashMap<Source, Double> sourceCosts = soln.getSourceCosts();
            for (Source src : sourceCaptureAmounts.keySet()) {
                bw.write(src.getLabel() + ",");
                bw.write(sourceCaptureAmounts.get(src) + ",");
                bw.write(sourceCosts.get(src) + "\n");
            }
            bw.write("\n");

            bw.write("Sink,Storage Amount (MTCO2/yr),Storage Cost ($M/yr)\n");
            HashMap<Sink, Double> sinkStorageAmounts = soln.getSinkStorageAmounts();
            HashMap<Sink, Double> sinkCosts = soln.getSinkCosts();
            for (Sink snk : sinkStorageAmounts.keySet()) {
                bw.write(snk.getLabel() + ",");
                bw.write(sinkStorageAmounts.get(snk) + ",");
                bw.write(sinkCosts.get(snk) + "\n");
            }

            bw.write("\n");

            bw.write("Edge Source,Edge Sink,Amount (MTCO2/yr),Transport Cost ($M/yr)\n");
            HashMap<Edge, Double> edgeTransportAmounts = soln.getEdgeTransportAmounts();
            HashMap<Edge, Double> edgeCosts = soln.getEdgeCosts();
            for (Edge edg : edgeTransportAmounts.keySet()) {
                bw.write(edg.v1 + "," + edg.v2 + ",");
                bw.write(edgeTransportAmounts.get(edg) + ",");
                bw.write(edgeCosts.get(edg) + "\n");
            }
        } catch (IOException e) {
            System.out.println(e.getMessage());
        }
    }

    public void makePriceAggregationFile(String path, String content) {
        try (BufferedWriter bw = new BufferedWriter(new FileWriter(path))) {
            bw.write(content);
        } catch (IOException e) {
            System.out.println(e.getMessage());
        }
    }

    public void makeGenerateFile(String path, Solution soln) {
        File newDir = new File(path + "/genFiles");
        if (true) {
            newDir.mkdir();
            Source[] sources = data.getSources();
            Sink[] sinks = data.getSinks();
            HashMap<Source, Double> sourceCaptureAmounts = soln.getSourceCaptureAmounts();
            HashMap<Sink, Double> sinkStorageAmounts = soln.getSinkStorageAmounts();
            HashMap<Edge, Double> edgeTransportAmounts = soln.getEdgeTransportAmounts();
            HashMap<Edge, int[]> graphEdgeRoutes = data.getGraphEdgeRoutes();

            // Make Sources.
            try (BufferedWriter bw = new BufferedWriter(new FileWriter(new File(newDir, "Sources.txt")))) {
                bw.write("ID,X,Y,CO2Cptrd,MxSpply,PieWdge,GensUsed,MaxGens,ActlCst,TtlCst,Name,Cell#\n");
                for (Source src : sources) {
                    bw.write(src.getLabel() + "," + data.cellToLatLon(src.getCellNum())[1] + "," + data.cellToLatLon(src.getCellNum())[0] + ",");
                    if (sourceCaptureAmounts.containsKey(src)) {
                        bw.write(sourceCaptureAmounts.get(src) + "," + src.getProductionRate() + "," + (src.getProductionRate() - sourceCaptureAmounts.get(src)));
                    } else {
                        bw.write("0," + src.getProductionRate() + "," + src.getProductionRate());
                    }
                    bw.write(",0,0,0,0,0,0\n");
                }
                bw.write("END");
            } catch (IOException e) {
                System.out.println(e.getMessage());
            }

            // Make Sinks.
            try (BufferedWriter bw = new BufferedWriter(new FileWriter(new File(newDir, "Sinks.txt")))) {
                bw.write("ID,X,Y,CO2Strd,MxStrg,PieWdge,WllsUsd,MxWlls,ActCst,TtlCst,Name,Cell#\n");
                for (Sink snk : sinks) {
                    bw.write(snk.getLabel() + "," + data.cellToLatLon(snk.getCellNum())[1] + "," + data.cellToLatLon(snk.getCellNum())[0] + ",");
                    if (sinkStorageAmounts.containsKey(snk)) {
                        bw.write(sinkStorageAmounts.get(snk) + "," + snk.getCapacity() + "," + (snk.getCapacity() - sinkStorageAmounts.get(snk)));
                    } else {
                        bw.write("0," + snk.getCapacity() + "," + snk.getCapacity());
                    }
                    bw.write(",0,0,0,0,0,0\n");
                }
                bw.write("END");
            } catch (IOException e) {
                System.out.println(e.getMessage());
            }

            // Make PipeDiameters.
            try (BufferedWriter bw = new BufferedWriter(new FileWriter(new File(newDir, "PipeDiameters.txt")))) {
                bw.write("ID,CapID,CapValue,Flow,Cost,LengKM,LengROW,LengCONS,Variable\n");
                for (Edge e : soln.getOpenedEdges()) {
                    bw.write("0,0,0," + edgeTransportAmounts.get(e) + ",0,0,0,0,0\n");
                    int[] route = graphEdgeRoutes.get(e);
                    for (int vertex : route) {
                        bw.write(round(data.cellToLatLon(vertex)[1], 5) + "," + round(data.cellToLatLon(vertex)[0], 5) + "\n");
                    }
                    bw.write("END\n");
                }
            } catch (IOException e) {
                System.out.println(e.getMessage());
            }
        }
    }

    // Download file from url
    public void downloadFile(String urlPath) {
        HttpURLConnection connection;

        try {
            URL url = new URL(urlPath);
            connection = (HttpURLConnection) url.openConnection();

            DateFormat dateFormat = new SimpleDateFormat("ddMMyyy-HHmmssss");
            Date date = new Date();
            String run = "run" + dateFormat.format(date);

            String directoryPath = basePath + "/" + dataset + "/Scenarios/" + scenario + "/Results/" + run;
            File directory = new File(directoryPath);
            directory.mkdir();

            // Copy MPS file into results file.
            String mipPath = basePath + "/" + dataset + "/Scenarios/" + scenario + "/MIP/mip.mps";
            Path from = Paths.get(mipPath);
            Path to = Paths.get(directoryPath + "/mip.mps");
            Files.copy(from, to, StandardCopyOption.REPLACE_EXISTING);

            FileOutputStream outputStream = null;
            InputStream inputStream = null;
            try {
                outputStream = new FileOutputStream(directoryPath + "/run0.sol");
                inputStream = connection.getInputStream();
                int BUFFER_SIZE = 10240;
                int bytesRead = -1;
                byte[] buffer = new byte[BUFFER_SIZE];
                while ((bytesRead = inputStream.read(buffer)) != -1) {
                    outputStream.write(buffer, 0, bytesRead);
                }
            } catch (FileNotFoundException ex) {
                ex.printStackTrace();
            } catch (IOException ex) {
                ex.printStackTrace();
            } finally {
                try {
                    if (inputStream != null) {
                        inputStream.close();
                    }
                } catch (IOException e) {
                }
                try {
                    if (outputStream != null) {
                        outputStream.close();
                    }
                } catch (IOException e) {
                }
            }
        } catch (MalformedURLException ex) {
            ex.printStackTrace();
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }
}
