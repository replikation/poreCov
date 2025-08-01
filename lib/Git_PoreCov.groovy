import java.net.URL
import groovy.json.JsonSlurper

class Git_PoreCov {
    /**
     * Gets the latest version of a GitHub repository
     * @return Latest version string
     */
    static String getLatestVersion() {
        try {
            final URL url = new URL("https://api.github.com/repos/replikation/poreCov/releases/latest")
            def response = url.getText(requestProperties: [
                Accept: 'application/json'
            ]);
            
            def json = new JsonSlurper().parseText(response);
            return json.tag_name ?: json.name;
        } catch (Exception e) {
            println "Error fetching version: ${e.message}";
            return 'Could not get version info';
        }
    }
}
