var wms_layers = [];


        var lyr_Positronnolabels_0 = new ol.layer.Tile({
            'title': 'Positron [no labels]',
            //'type': 'base',
            'opacity': 1.000000,
            
            
            source: new ol.source.XYZ({
    attributions: ' &middot; <a href="https://cartodb.com/basemaps/">Map tiles by CartoDB, under CC BY 3.0. Data by OpenStreetMap, under ODbL.</a>',
                url: 'https://a.basemaps.cartocdn.com/light_nolabels/{z}/{x}/{y}.png'
            })
        });
var format_Countries_1 = new ol.format.GeoJSON();
var features_Countries_1 = format_Countries_1.readFeatures(json_Countries_1, 
            {dataProjection: 'EPSG:4326', featureProjection: 'EPSG:3857'});
var jsonSource_Countries_1 = new ol.source.Vector({
    attributions: ' ',
});
jsonSource_Countries_1.addFeatures(features_Countries_1);
var lyr_Countries_1 = new ol.layer.Vector({
                declutter: false,
                source:jsonSource_Countries_1, 
                style: style_Countries_1,
                popuplayertitle: "Countries",
                interactive: true,
    title: 'Countries<br />\
    <img src="styles/legend/Countries_1_0.png" /> Cape Verde<br />\
    <img src="styles/legend/Countries_1_1.png" /> Cook Islands<br />\
    <img src="styles/legend/Countries_1_2.png" /> Cote d\'lvoire<br />\
    <img src="styles/legend/Countries_1_3.png" /> Gabon<br />\
    <img src="styles/legend/Countries_1_4.png" /> Greenland<br />\
    <img src="styles/legend/Countries_1_5.png" /> Guinea-Bissau <br />\
    <img src="styles/legend/Countries_1_6.png" /> Kiribati<br />\
    <img src="styles/legend/Countries_1_7.png" /> Madagascar<br />\
    <img src="styles/legend/Countries_1_8.png" /> Mauritania<br />\
    <img src="styles/legend/Countries_1_9.png" /> Mauritius<br />\
    <img src="styles/legend/Countries_1_10.png" /> São Tomé and Principe<br />\
    <img src="styles/legend/Countries_1_11.png" /> Senegal <br />\
    <img src="styles/legend/Countries_1_12.png" /> Seychelles<br />\
    <img src="styles/legend/Countries_1_13.png" /> The Gambia<br />'
        });

lyr_Positronnolabels_0.setVisible(true);lyr_Countries_1.setVisible(true);
var layersList = [lyr_Positronnolabels_0,lyr_Countries_1];
lyr_Countries_1.set('fieldAliases', {'COUNTRY': 'COUNTRY', 'ISO': 'ISO', 'COUNTRYAFF': 'COUNTRYAFF', 'AFF_ISO': 'AFF_ISO', 'Name': 'Name', 'Area of EEZ (km^2)': 'Area of EEZ (km^2)', 'Population Density (people per km^2 of land area) 2021': 'Population Density (people per km^2 of land area) 2021', 'Area of Land Mass (excluding inland water) (km^2 )': 'Area of Land Mass (excluding inland water) (km^2 )', 'Colonial Power': 'Colonial Power', 'Year of Independence': 'Year of Independence', 'Size of National Fleet (vessels) 2023': 'Size of National Fleet (vessels) 2023', '3-Letter Country Code': '3-Letter Country Code', 'Alpha-3 Code': 'Alpha-3 Code', });
lyr_Countries_1.set('fieldImages', {'COUNTRY': 'TextEdit', 'ISO': 'TextEdit', 'COUNTRYAFF': 'TextEdit', 'AFF_ISO': 'TextEdit', 'Name': 'TextEdit', 'Area of EEZ (km^2)': 'Range', 'Population Density (people per km^2 of land area) 2021': 'TextEdit', 'Area of Land Mass (excluding inland water) (km^2 )': 'Range', 'Colonial Power': 'TextEdit', 'Year of Independence': 'DateTime', 'Size of National Fleet (vessels) 2023': 'Range', '3-Letter Country Code': 'Range', 'Alpha-3 Code': 'TextEdit', });
lyr_Countries_1.set('fieldLabels', {'COUNTRY': 'hidden field', 'ISO': 'hidden field', 'COUNTRYAFF': 'hidden field', 'AFF_ISO': 'hidden field', 'Name': 'no label', 'Area of EEZ (km^2)': 'inline label - visible with data', 'Population Density (people per km^2 of land area) 2021': 'hidden field', 'Area of Land Mass (excluding inland water) (km^2 )': 'hidden field', 'Colonial Power': 'inline label - visible with data', 'Year of Independence': 'inline label - visible with data', 'Size of National Fleet (vessels) 2023': 'inline label - visible with data', '3-Letter Country Code': 'hidden field', 'Alpha-3 Code': 'hidden field', });
lyr_Countries_1.on('precompose', function(evt) {
    evt.context.globalCompositeOperation = 'normal';
});