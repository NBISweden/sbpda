import static qupath.lib.gui.scripting.QPEx.*
import javafx.application.Platform
import qupath.lib.images.ImageData
import qupath.lib.images.servers.ImageServerProvider
import qupath.lib.images.servers.TransformedServerBuilder


// Please provide path to output directory where split images for each channel
PATH_TO_OUTPUT_DIR = '/Users/Desktop/'

// Get the current image
def imageData = getCurrentImageData()
def server = imageData.getServer()
describe(server)
// Define the output directory
def outputDir = buildFilePath(PATH_TO_OUTPUT_DIR, 'exported_channels')
mkdirs(outputDir)
// Loop through each channel
def channelNames = getCurrentServer().getMetadata().getChannels().collect { c -> c.name }
print('There are ' + channelNames.size() + ' channels/markers in this image')



// Loop through each channel
channelNames.eachWithIndex { channelName, index ->
    // Print the channel name
    def outputPath = buildFilePath(outputDir, "${channelName}.tif")
    print("Channel " + index + ": " + channelName + ' ==> '+ outputPath)
    server2 = new TransformedServerBuilder(server)
    .extractChannels(index)
    .build()
    
//    def outputPath = buildFilePath(outputDir, "extracted_channel.tif")
    def requestFull = RegionRequest.createInstance(server2, 1)
    writeImageRegion(server2, requestFull, outputPath)
}

print("Image saved")