package mvieths.genomics;

import java.util.ArrayList;

import javafx.application.Application;
import javafx.scene.Scene;
import javafx.scene.chart.LineChart;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.XYChart;
import javafx.stage.Stage;

public class GCGeneChart extends Application {
    static ArrayList<GCATWindow> windows;
    static ArrayList<GeneWindow> genes;
    static Chromosome            chr;
    static long                  window;

    @Override
    public void start(Stage stage) {
        stage.setTitle("Nucleotide Percentages in Chromosome 13");
        final NumberAxis xAxis = new NumberAxis(46000000, 47000000, 100000);
        final NumberAxis yAxis = new NumberAxis(0, 75.0, 5.0);
        xAxis.setLabel("Nucleotide Position (Window Size = " + window + ")");
        yAxis.setLabel("Percentage");

        final LineChart<Number, Number> lineChart = new LineChart<Number, Number>(xAxis, yAxis);
        lineChart.setTitle("Nucleotide Percentages in Chromosome 13");

        XYChart.Series<Number, Number> gcSeries = new XYChart.Series<Number, Number>();
        gcSeries.setName("GC Content");

        for (GCATWindow window : windows) {
            gcSeries.getData().add(
                    new XYChart.Data<Number, Number>(window.getWindowNumber(), (long) window.getGCPercent()));
        }

        lineChart.getData().add(gcSeries);

        for (int i = 0; i < genes.size(); i++) {
            GeneWindow gene = genes.get(i);
            double yCoord = 50.0 + (i % 25);

            XYChart.Series<Number, Number> geneSeries = new XYChart.Series<Number, Number>();
            geneSeries.setName(gene.getName());
            geneSeries.getData().add(new XYChart.Data<Number, Number>(gene.getStart(), yCoord));
            geneSeries.getData().add(new XYChart.Data<Number, Number>(gene.getEnd(), yCoord));
            lineChart.getData().add(geneSeries);
        }

        Scene scene = new Scene(lineChart, 1200, 800);

        stage.setScene(scene);
        stage.show();
    }

    public static void main(String[] args) {
        chr = new Chromosome(args[0], args[1]);
        chr.calcContent();
        chr.calcGenes();
        long chrLength = chr.getChromosomeLength();
        // We want a readable chart, so limit it to 50 data points
        window = chrLength / 1000;
        chr.calcContentByWindow(window);

        windows = chr.getWindows();
        genes = chr.getGenes();

        launch(args);
    }
}