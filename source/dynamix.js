//Dynamix Desktop Environment
var dynamix_version="1.0";
var swidth=width;
var sheight=height;
var package_names=[]; //Good Packages (Version Keys)
var window_open=false;
var rmxbar=false; //Remove Multi-use X Bar
var opensans=createFont("Google Sans"); //Google Sans Font (Default)
var monospace=createFont("monospace"); //Monospace Font
var lato=createFont("Lato"); //Lato Font
var roboto=createFont("Roboto"); //Roboto Font
var times_new_roman=createFont("Times New Roman"); //Times New Roman Font
var ubuntu=createFont("Ubuntu"); //Ubuntu Font
var arial=createFont("Arial"); //Arial Font
var georgia=createFont("Georgia"); //Georgia Font
var clean_dsk=function(){
frameRate(60);
smooth();
lights();
noStroke();
}; //Cleans up shapes, etc
var window=function(txt,r,g,b,a){
    if(window_open===true){
fill(r,g,b,a);
rect(0,sheight/15,swidth+3,sheight-swidth/1000);
if(rmxbar===false){
fill(0, 0, 0,80);
rect(swidth-200,sheight/15.2,swidth/2.5,sheight/15);
fill(255, 255, 255);
textFont(opensans,swidth/17);
text("X",swidth-25,sheight/8.5);
textFont(opensans,swidth/25);
text(txt,swidth-200,sheight/8.5);
if(mouseX>=swidth-34&&mouseY>=sheight/10&&mouseY<=sheight/10&&mousePressed){
    window_open=false;
}
} else{}
 } else{}
}; //Opens a new window (title,color of window)
var scriptspider=function(){
    println("---ScriptSpider Build 0000A---\n");
    println("Checking List...");
    println("There are "+package_names.length+" Good Packages to Search.");
    //Sorry No Apps Yet :'(
    
}; //Finds a package from a list of known goods (Apps not libaries)
var cursorset=function(num){
if(num===0){
cursor("NONE");
 } //No Cursor
 if(num===1){
 cursor("NONE");
 stroke(0, 0, 0);
 fill(245, 245, 245);
 ellipse(mouseX,mouseY,20,20);
 clean_dsk();
 } //Default Cursor
 if(num===1.5){
 cursor("NONE");
 stroke(0, 0, 0);
 fill(245, 245, 245);
 ellipse(mouseX,mouseY,20,20);
 clean_dsk();
 fill(255, 0, 0);
 text(round(swidth/mouseX)+" , "+round(sheight/mouseY),mouseX-10,mouseY-10);
 } //Debug Cursor (For Development)
}; //Cursor Setup
var dynamix=function(){
    scriptspider();
draw= function() {
    frameRate(60);
    var m=minute();
    var h=hour();
    var dm,dh; //Display Clock Settings
    if(m<10){
        dm="0"+m;
    } 
    else{dm=m;
    
    }
    if(h<10){
        dh="0"+hour;
    }
    else{
        dh=h;
    }
    //Basic UI Stuff
 background(199, 199, 199);
 clean_dsk(); //Anti-Alising Stuff
 fill(102, 102, 102);
 rect(-5,sheight/10000,swidth+8,sheight/15);
 fill(199, 199, 199);
 textFont(opensans,sheight/20);
 text(dh+":"+dm,swidth/1.15,sheight/20);
    
 
    cursorset(1); //Always Goes Last <----
 };
}; //Main Dynamix Desktop

