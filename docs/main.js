
var info, detail, current_protein, json_data;
var DATA_FILE = "./challenge.json";
var WIDTH = 1000;
var select, detail, checkboxes;

function setup() {
    info = $('#info');
    detail = $('#protein_detail');
    current_protein = null;
    info.html("Javascript loaded.");
    $.getJSON(DATA_FILE, process_json).fail(on_load_failure);
};

var on_load_failure = function() {
    alert("Could not load local JSON data.\n" +
            "You may need to run a web server to avoid cross origin restrictions.")
};

var process_json = function(data) {
    debugger;
    json_data = data;
    var p = "loaded proteins: ";
    var selector = $("#protein_select");
    detail = $('#protein_detail');
    selector.empty()
    $("<span>Protein: </span>").appendTo(selector);
    select = $('<select id="protein_selection"/>').appendTo(selector);
    for (var name in data) {
        var option = $(`<option value="${name}">${name}</option>`).appendTo(select);
        p += name + "; ";
        if (!current_protein) {
            current_protein = data[name]
        }
    }
    select.selectmenu();
    select.on('selectmenuchange', draw_selection);
    //$("#protein_selection").on('selectmenuchange', draw_selection);
    info.html(p);
    draw_selection();
};

var other_classification = "(other)"

var draw_selection = function() {
    debugger;
    var name = select.val();
    current_protein = json_data[name];
    detail.empty();
    checkboxes = [];
    var classifications = current_protein.classifications;
    for (var i=0; i<classifications.length; i++) {
        var classification = classifications[i].trim();
        if (classification.length < 1) {
            classification = other_classification;
        }
        var cbdiv = $("<div/>").appendTo(detail);
        var cb = $(`<input type="checkbox" value="${classification}" checked/>`).appendTo(cbdiv);
        $(`<span> ${classification} </span>`).appendTo(cbdiv);
        checkboxes.push(cb);
    }
    prepare_canvas();
    draw_protein();
    //nd_frame.fit(0.8)
    nd_frame.orbit_all(current_protein.radius, current_protein.center);
    debugger;
    nd_frame.rotate_shift(current_protein.center, current_protein.radius, [5,0]);
};

var target, frame, nd_frame;

var prepare_canvas = function() {
    target = $("<div/>").appendTo(detail);
    var canvas_config = {
        width: WIDTH,
        height: WIDTH,
    };
    var r2 = current_protein.radius * 0.4;
    target.dual_canvas_helper(canvas_config);
    var [cx, cy, cz] = current_protein.center;
    frame = target.frame_region(
        0, 0, WIDTH, WIDTH,
        cx-r2, cy-r2, cx+r2, cy+r2);
    nd_frame = target.nd_frame({
        dedicated_frame: frame,
    });
};

var draw_protein = function() {
    annotation_radius = 10;
    var locations = current_protein.locations;
    //var transparent = "rgba(0,0,0,0)";
    nd_frame.polygon({
        locations: locations,
        close: false,
        fill: false,
    });
    var residues = current_protein.residues;
    for (var rname in residues) {
        var r = residues[rname];
        var res = r.RES;
        var acid = Amino_acids[res];
        var [R,G,B] = acid.color;
        var color = "rgb(" + R + "," + G + "," + B + ")";
        var semi_transparent = "rgba(" + R + "," + G + "," + B + ", 0.5)";
        nd_frame.polygon({
            locations: [r.C.pos, r.N.pos, r.CA.pos],
            color: color,
            fill: false,
            lineWidth: 2,
        });
        if (r.CB) {
            nd_frame.polygon({
                locations: [r.CA.pos, r.CB.pos],
                close: false,
                fill: false,
                color: color,
                lineWidth: 3,
            });
        }
        var annotations = r.annotations;
        var n_ann = annotations.length;
        if (n_ann > 0) {
            // add an annotation marker
            nd_frame.circle({
                location:r.C.pos, 
                r: annotation_radius,
                color:color,
                fill:false,
                lineWidth:2,
            });
            var rcircle = nd_frame.circle({
                location:r.C.pos, 
                r: annotation_radius,
                color:semi_transparent,
                name: "residue_" + rname,
            });
            rcircle.on("mouseover", rcircle_mouseover)
        }
    }
};

var rcircle_mouseover = function(event) {
    debugger;
    var circle_name = event.canvas_name;
    var residue_name = circle_name.split("_")[1]
    info.html("Hover over residue: " + residue_name)
}

Amino_acids = {
'A': {'name': 'Alanine',
  'abbrev3': 'Ala',
  'abbrev1': 'A',
  'color': [140, 213, 140]},
 'R': {'name': 'Arginine',
  'abbrev3': 'Arg',
  'abbrev1': 'R',
  'color': [0, 0, 124]},
 'N': {'name': 'Asparagine',
  'abbrev3': 'Asn',
  'abbrev1': 'N',
  'color': [213, 124, 112]},
 'D': {'name': 'Aspartic acid',
  'abbrev3': 'Asp',
  'abbrev1': 'D',
  'color': [160, 0, 66]},
 'C': {'name': 'Cysteine',
  'abbrev3': 'Cys',
  'abbrev1': 'C',
  'color': [191, 191, 91]},
 'E': {'name': 'Glutamic acid',
  'abbrev3': 'Glu',
  'abbrev1': 'E',
  'color': [102, 0, 0]},
 'Q': {'name': 'Glutamine',
  'abbrev3': 'Gln',
  'abbrev1': 'Q',
  'color': [213, 76, 76]},
 'G': {'name': 'Glycine',
  'abbrev3': 'Gly',
  'abbrev1': 'G',
  'color': [187, 187, 187]},
 'H': {'name': 'Histidine',
  'abbrev3': 'His',
  'abbrev1': 'H',
  'color': [112, 112, 213]},
 'I': {'name': 'Isoleucine',
  'abbrev3': 'Ile',
  'abbrev1': 'I',
  'color': [0, 76, 0]},
 'L': {'name': 'Leucine',
  'abbrev3': 'Leu',
  'abbrev1': 'L',
  'color': [69, 94, 69]},
 'K': {'name': 'Lysine',
  'abbrev3': 'Lys',
  'abbrev1': 'K',
  'color': [71, 71, 184]},
 'M': {'name': 'Methionine',
  'abbrev3': 'Met',
  'abbrev1': 'M',
  'color': [184, 160, 66]},
 'F': {'name': 'Phenylalanine',
  'abbrev3': 'Phe',
  'abbrev1': 'F',
  'color': [83, 76, 66]},
 'P': {'name': 'Proline',
  'abbrev3': 'Pro',
  'abbrev1': 'P',
  'color': [82, 82, 82]},
 'S': {'name': 'Serine',
  'abbrev3': 'Ser',
  'abbrev1': 'S',
  'color': [213, 112, 66]},
 'T': {'name': 'Threonine',
  'abbrev3': 'Thr',
  'abbrev1': 'T',
  'color': [184, 76, 0]},
 'W': {'name': 'Tryptophan',
  'abbrev3': 'Trp',
  'abbrev1': 'W',
  'color': [79, 70, 0]},
 'Y': {'name': 'Tyrosine',
  'abbrev3': 'Tyr',
  'abbrev1': 'Y',
  'color': [140, 112, 76]},
 'V': {'name': 'Valine',
  'abbrev3': 'Val',
  'abbrev1': 'V',
  'color': [213, 140, 213]}}