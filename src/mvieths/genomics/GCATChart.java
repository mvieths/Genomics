package mvieths.genomics;

import java.util.ArrayList;

import javafx.application.Application;
import javafx.scene.Scene;
import javafx.scene.chart.CategoryAxis;
import javafx.scene.chart.LineChart;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.XYChart;
import javafx.stage.Stage;

public class GCATChart extends Application {
    static ArrayList<GCATWindow> windows;
    static Chromosome            chr;
    static long                  window;

    @Override
    public void start(Stage stage) {
        stage.setTitle("G/C/A/T Percentages in Chromosome 13");
        final CategoryAxis xAxis = new CategoryAxis();
        final NumberAxis yAxis = new NumberAxis(0, 50.0, 5.0);
        xAxis.setLabel("Nucleotide Position (Window Size = " + window + ")");
        yAxis.setLabel("Percentage");

        final LineChart<String, Number> lineChart = new LineChart<String, Number>(xAxis, yAxis);
        lineChart.setTitle("G/C/A/T Percentages in Chromosome 13");
        XYChart.Series<String, Number> gSeries = new XYChart.Series<String, Number>();
        gSeries.setName("G Content");
        XYChart.Series<String, Number> cSeries = new XYChart.Series<String, Number>();
        cSeries.setName("C Content");
        XYChart.Series<String, Number> aSeries = new XYChart.Series<String, Number>();
        aSeries.setName("A Content");
        XYChart.Series<String, Number> tSeries = new XYChart.Series<String, Number>();
        tSeries.setName("T Content");

        for (GCATWindow window : windows) {
            gSeries.getData().add(new XYChart.Data<String, Number>(window.getWindowNumber(), window.getgPercent()));
            cSeries.getData().add(new XYChart.Data<String, Number>(window.getWindowNumber(), window.getcPercent()));
            aSeries.getData().add(new XYChart.Data<String, Number>(window.getWindowNumber(), window.getaPercent()));
            tSeries.getData().add(new XYChart.Data<String, Number>(window.getWindowNumber(), window.gettPercent()));
        }

        Scene scene = new Scene(lineChart, 1200, 800);

        lineChart.getData().addAll(gSeries, cSeries, aSeries, tSeries);

        stage.setScene(scene);
        stage.show();
    }

    public static void main(String[] args) {
        chr = new Chromosome("C:\\Users\\Foeclan\\Dropbox\\Genomics\\chr13.fa");
        chr.calcContent();
        chr.printGCPercent();
        chr.printAllPercents();
        long chrLength = chr.getChromosomeLength();
        // We want a readable chart, so limit it to 50 data points
        window = chrLength / 50;
        chr.calcContentByWindow(window);

        windows = chr.getWindows();

        launch(args);
    }
}