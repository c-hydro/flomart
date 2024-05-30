<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>
<qgis version="3.4.13-Madeira" maxScale="0" minScale="1e+08" hasScaleBasedVisibilityFlag="0" styleCategories="AllStyleCategories">
  <flags>
    <Identifiable>1</Identifiable>
    <Removable>1</Removable>
    <Searchable>1</Searchable>
  </flags>
  <customproperties>
    <property key="WMSBackgroundLayer" value="false"/>
    <property key="WMSPublishDataSourceUrl" value="false"/>
    <property key="embeddedWidgets/count" value="0"/>
    <property key="identify/format" value="Value"/>
  </customproperties>
  <pipe>
    <rasterrenderer alphaBand="-1" opacity="0.655" type="singlebandpseudocolor" classificationMin="0" classificationMax="1" band="1">
      <rasterTransparency/>
      <minMaxOrigin>
        <limits>None</limits>
        <extent>WholeRaster</extent>
        <statAccuracy>Estimated</statAccuracy>
        <cumulativeCutLower>0.02</cumulativeCutLower>
        <cumulativeCutUpper>0.98</cumulativeCutUpper>
        <stdDevFactor>2</stdDevFactor>
      </minMaxOrigin>
      <rastershader>
        <colorrampshader colorRampType="INTERPOLATED" clip="0" classificationMode="1">
          <colorramp name="[source]" type="gradient">
            <prop v="247,251,255,255" k="color1"/>
            <prop v="8,48,107,255" k="color2"/>
            <prop v="0" k="discrete"/>
            <prop v="gradient" k="rampType"/>
            <prop v="0.13;222,235,247,255:0.26;198,219,239,255:0.39;158,202,225,255:0.52;107,174,214,255:0.65;66,146,198,255:0.78;33,113,181,255:0.9;8,81,156,255" k="stops"/>
          </colorramp>
          <item alpha="255" color="#f7fbff" label="0" value="0"/>
          <item alpha="255" color="#deebf7" label="0.13" value="0.13"/>
          <item alpha="255" color="#c6dbef" label="0.26" value="0.26"/>
          <item alpha="255" color="#9ecae1" label="0.39" value="0.39"/>
          <item alpha="255" color="#6baed6" label="0.52" value="0.52"/>
          <item alpha="255" color="#4292c6" label="0.65" value="0.65"/>
          <item alpha="255" color="#2171b5" label="0.78" value="0.78"/>
          <item alpha="255" color="#08519c" label="0.9" value="0.9"/>
          <item alpha="255" color="#08306b" label="1" value="1"/>
        </colorrampshader>
      </rastershader>
    </rasterrenderer>
    <brightnesscontrast brightness="0" contrast="0"/>
    <huesaturation grayscaleMode="0" colorizeBlue="128" saturation="0" colorizeRed="255" colorizeGreen="128" colorizeOn="0" colorizeStrength="100"/>
    <rasterresampler maxOversampling="2"/>
  </pipe>
  <blendMode>0</blendMode>
</qgis>
