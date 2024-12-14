proc make_frame { filename id_num rotation} {
  mol new $filename type {mol2}
  mol modstyle 0 $id_num CPK 2.00000 0.000000 22.000000 15.000000
  mol modcolor 0 $id_num Charge
  mol scaleminmax $id_num 0 -1.00000 1.00000
  mol modmaterial 0 $id_num Glossy

  axes location Off
  #menu colorscalebar on
  #rotate x by -10.5
  #rotate x by 10.5
  #rotate y by $rotation


  set outputfile [file rootname $filename].bmp
  render TachyonInternal $outputfile
  mol delete $id_num
}

proc main_loop { path } {
  puts "start"
  set mol2_list [glob $path/*.mol2]
  puts "after mol2_list"
  puts $mol2_list
  set sorted_list [lsort -command {apply {{a b} {
    # Extract the number from each filename by stripping the path and extension
    if {![regexp {(\d+)\.mol2$} [file tail $a] _ a_num] || ![regexp {(\d+)\.mol2$} [file tail $b] _ b_num]} {
        # If either filename doesn't match the expected pattern, treat as equal
        return 0
    }

    # Return the comparison result
    expr {$a_num - $b_num}
}}} $mol2_list]
  puts "after sort list"

  set count 0
  set rotate_amount 0

  foreach mol2 $sorted_list { 
    make_frame $mol2 $count $rotate_amount
    puts "Current value of count: $count"
    incr count
    #incr rotate_amount 3

  }
} 