import dataStore.DataStorer;
import solver.Solver;

public class Main {

    public static void main(String[] args) {
        
        String basepath = "./Datasets";
        String dataset = "Lower48US";
        String scenario = "scenario1";
        DataStorer data = new DataStorer(basepath, dataset, scenario);
        Solver solver = new Solver(data);
        data.setSolver(solver);
        data.loadNetworkCosts();
    }
}
