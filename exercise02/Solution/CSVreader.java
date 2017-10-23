import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

public class CSVreader {
	private static String separator = ";";

	public static List<List<Double>> readRecords(Path path) {

		List<List<Double>> values;

		try (BufferedReader br = Files.newBufferedReader(path)) {
			values = br.lines().map(line -> Arrays.asList(line.split(separator)).stream().map(Double::parseDouble)
					.collect(Collectors.toList())).collect(Collectors.toList());

			return values;
		} catch (IOException e) {
			e.printStackTrace();
			return null;
		}
	}

	private static List<Double> mean(List<List<Double>> values) {
		List<Double> avg = new ArrayList<Double>();

		for (int i = 0; i < values.get(1).size(); i++) {
			Double tmpAvg = 0.;
			for (int j = 0; j < values.size(); j++) {
				tmpAvg += values.get(j).get(i);
			}
			tmpAvg /= values.size();
			avg.add(tmpAvg);
		}
		return avg;
	}

	public static List<Double> stdDeviation(List<List<Double>> values, List<Double> mean) {
		List<Double> stdDeviation = new ArrayList<Double>();

		for (int i = 0; i < values.get(1).size(); i++) {
			Double numj = 0.;
			Double num = 0.;
			Double meanTmp = mean.get(i);
			for (int j = 0; j < values.size(); j++) {
				numj = Math.pow(values.get(j).get(i) - meanTmp, 2);
				num += numj;
			}

			stdDeviation.add(Math.sqrt(num / values.size()));
		}

		return stdDeviation;

	}

	public static String listToCSV(List<Double> list) {
		String out = "";

		for (Double tmp : list) {
			out += tmp + separator;
		}

		return out.substring(0, out.length() - 1);
	}

	public static void main(String[] args) {

		File workingDir = new File("/Users/Philipp/ParSys_LP/exercise02/Solution/part 1/csv");

		for (File file : workingDir.listFiles()) {

			if (file.isHidden() || file.getName().equals("avgResults.csv") || file.getName().equals("stdDevResults.csv"))
				continue;

			calculateForFile(file.toPath());
		}

	}

	private static void calculateForFile(Path input) {
		List<List<Double>> values = readRecords(input);
		// System.out.println(values);

		List<Double> mean = mean(values);
		List<Double> stdDev = stdDeviation(values, mean);

		String fileName = input.getFileName().toString().substring(0, input.getFileName().toString().length() - 4);
		System.out.println(fileName + "Avg;" + listToCSV(mean));
		System.out.println(fileName + "StdDev;" + listToCSV(stdDev));
		
		try {
			Files.write(Paths.get(input.getParent() +"/avgResults.csv"),
					(fileName + ";" + listToCSV(mean) + "\n").getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE,
					StandardOpenOption.APPEND);
			
			Files.write(Paths.get(input.getParent() +"/stdDevResults.csv"),
					(fileName + ";" + listToCSV(stdDev) + "\n").getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE,
					StandardOpenOption.APPEND);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

}
