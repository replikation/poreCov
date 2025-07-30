import java.net.URL
import java.net.URLConnection

class Hgdownload_lcs {
    static boolean IsLcsAvailable() {
        try {
            final URL url = new URL("https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/public-latest.version.txt");
            final URLConnection conn = url.openConnection();
            conn.connect();
            conn.getInputStream().close();
            return true;
        } catch (MalformedURLException e) {
            return false;
        } catch (IOException e) {
            return false;
        }
    }
}