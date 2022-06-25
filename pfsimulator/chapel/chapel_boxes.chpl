use CTypes;
extern type Point = 3*int(32);

extern record Box {
    var lo: Point;
    var up: Point;
};

extern record BoxArray {
  var size: int;
  var boxes: [0..<size] Box;
};

export proc print_boxes_from_surface_boxes(ref surface_boxes: [0..<8] BoxArray, PV_f: int) {
  writeln("Printing boxes from the surface box list.");
  var curr_array: BoxArray = surface_boxes[0];

  //var boxes = curr_array.boxes;
  //for box in boxes {
  //  print_box(box);
  //}
  return 0;
}



export proc print_box(in b: Box) {
    writeln("box: ", b.lo, " ", b.up);
}

export proc chapel_print() {
    writeln("Printing from chapel code.");
}

export proc print_boxes(ref boxes: [] Box) {
    for box in boxes {
      print_box(box);
    }
}
export proc print_boxes_from_box_array(ref box_array: BoxArray) {
  var box = box_array.boxes;
} 
proc GrGeomSurfLoopBoxes_phase_rel_perm(ix: int, iy: int, iz: int, 
                                   nx: int, ny: int, nz: int, 
                                   fdir: [] int, ref boxes: [] Box) {

}

