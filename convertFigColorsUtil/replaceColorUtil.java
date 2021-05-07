//
// replacing colors for an image according to a user defined color replace table
//
// usage (args):
//
// 1) "getColors in.png [considerGrayscale]" -- get initial colors
//
// 2) "replaceColors in.png out.png colorMap.txt" -- replace colors; colorMap.txt looks like:
//
//   255,   0, 255 -> 210, 210, 210
//   255,   0,   0 -> 150, 150, 150
//     0,   0, 255 ->  90,  90,  90
//


import java.awt.Color;
import java.awt.image.BufferedImage;

import java.io.File;
import java.io.IOException;

import java.nio.file.Files;
import java.nio.file.Paths;

import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Stream;

import javax.imageio.ImageIO;



public class replaceColorUtil {

    private final static String LS = System.lineSeparator();

    private BufferedImage img = null;
    private int w = 0, h = 0;

    Set<Integer> rgbSet;

    private void load(String path) throws IOException {

        img = ImageIO.read(new File(path));
        w = img.getWidth();
        h = img.getHeight();
    }

    private void save(String path) throws IOException {

        ImageIO.write(img, "png", new File(path));
    }

    public void colorList(boolean considerGrayscale) {

        rgbSet = new TreeSet<>();

        for (int y = 0; y < h; ++y) {
            for (int x = 0; x < w; ++x) {
                rgbSet.add(img.getRGB(x, y));
            }
        }

        for (int c: rgbSet) {

            Color col = new Color(c);
            int r = col.getRed(), g = col.getGreen(), b = col.getBlue();

            boolean grayscale = ((r == g) && (g == b));

            if (considerGrayscale || !grayscale) {
                System.out.println(r + " " + g + " " + b);
            }
        }
    }

    private void replaceColor(int from, int to) {

        for (int y = 0; y < h; ++y) {
            for (int x = 0; x < w; ++x) {
                if (from == img.getRGB(x, y)) {
                    img.setRGB(x, y, to);
                }
            }
        }
    }

    private int[] parseColors(String mapLine) {

        String formatErr =
                "error: invalid map line format, s = \"" + mapLine + "\"" +
                LS + "expected: \"<R1>, <G1>, <B1> -> <R2>, <G2>, <B2>\";" +
                " the values should be in [0..255] range" + LS +
                "skipping this line...";

        String tmp[] = mapLine.trim().split("->");
        if (tmp.length != 2) {
            System.err.println(formatErr);
            return null;
        }

        int cols[] = new int[2];

        for (int i = 0; i < 2; ++i) {

            String c[] = tmp[i].trim().split(",");
            if (c.length != 3) {
                System.err.println(formatErr);
                return null;
            }

            int rgb[] = new int[3];
            for (int j = 0; j < 3; ++j) {

                int cc = Integer.parseInt(c[j].trim());
                if ((cc < 0) || (cc > 255)) {
                    System.err.println(formatErr);
                    return null;
                }
                rgb[j] = cc;
            }

            cols[i] = new Color(rgb[0], rgb[1], rgb[2]).getRGB();
        }

        return cols;
    }

    private void applyMapLine(String mapLine) {

        String s = mapLine.trim();
        if (s.isEmpty() || s.startsWith("#")) { return; }
        int cols[] = parseColors(s);
        replaceColor(cols[0], cols[1]);
    }



    public static void main(String args[]) throws Exception {

        String usage =
                "expected args: " + LS +
                "  getColors in.png [considerGrayscale]" + LS + "or" + LS +
                "  replaceColors in.png out.png colorMap.txt";

        if (args.length < 1) {
            System.err.println(usage);
            return;
        }

        replaceColorUtil rc = new replaceColorUtil();

        String cmd = args[0];

        if (cmd.equals("getColors")) {

            if (args.length < 2) {
                System.err.println("please provide input file path");
                return;
            }

            boolean considerGrayscale = (args.length > 2);

            rc.load(args[1]);
            rc.colorList(considerGrayscale);

        } else if (cmd.equals("replaceColors")) {

            if (args.length < 4) {
                System.err.println("expected replaceColors args: " +
                        "inFile.png outFile.png colorMap.txt");
                return;
            }

            String in = args[1], out = args[2], map = args[3];

            rc.load(in);
            try (Stream<String> lines = Files.lines(Paths.get(map))) {
                lines.forEach(rc::applyMapLine);
            }
            rc.save(out);

        } else {

            System.err.println("unknown command: " + cmd);
        }
    }
}
