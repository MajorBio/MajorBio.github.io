# 联系我们 


<!--link href="https://maxcdn.bootstrapcdn.com/font-awesome/4.7.0/css/font-awesome.min.css" rel="stylesheet" integrity="sha384-wvfXpqpZZVQGK6TAh5PVlGOfQNHSoD2xbE+QkPxCAFlNEevoEH3Sl0sibVcOQVnN" crossorigin="anonymous">
<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/animate.css/3.5.2/animate.css" /-->
<link rel="stylesheet" href="styles.css">

<br>


<div class="contact">
  <form action="https://formspree.io/meng.luo@majorbio.com" method="POST" enctype="multipart/form-data" name="EmailForm">
    <p>
      <label>Name *</label>
      <input type="text" name="Name" required>
    </p>
    <p>
      <label>Current University Affiliation *</label>
      <input type="text" name="Current University Affiliation" required>
    </p>
      <label>Current Academic Standing *</label>
    <p>
    <div class="custom-select" style="width:100%;" name="Current Academic Standing">
      <select id="Current Academic Standing" name="Current Academic Standing">
        <option value="Undergraduate">Undergraduate</option>
        <option value="Undergraduate">Undergraduate</option>
        <option value="Graduate">Graduate</option>
        <option value="PhD">PhD</option>
        <option value="Post-Doctoral Candidate">Post-Doctoral Candidate</option>
        <option value="Other" >Other</option>
      </select>
    </div>
    </p>
    <p>
      <label>Email Address *</label>
      <input type="email" name="Email Address" required>
    </p>
    <p>
      <label>Phone Number *</label>
      <input type="text" name="Phone Number" required>
    </p>
    <p class="full">
      <label>Brief Biography *</label>
      <textarea style="width:83%;" name="Brief Biography" rows="5" required></textarea>
    </p>
    <p class="full">
    <button style="width:90%;" align="center" >Submit</button>
    </p>
  </form>
</div>


<script>
var x, i, j, l, ll, selElmnt, a, b, c;
/*look for any elements with the class "custom-select":*/
x = document.getElementsByClassName("custom-select");
l = x.length;
for (i = 0; i < l; i++) {
  selElmnt = x[i].getElementsByTagName("select")[0];
  ll = selElmnt.length;
  /*for each element, create a new DIV that will act as the selected item:*/
  a = document.createElement("DIV");
  a.setAttribute("class", "select-selected");
  a.innerHTML = selElmnt.options[selElmnt.selectedIndex].innerHTML;
  x[i].appendChild(a);
  /*for each element, create a new DIV that will contain the option list:*/
  b = document.createElement("DIV");
  b.setAttribute("class", "select-items select-hide");
  for (j = 1; j < ll; j++) {
    /*for each option in the original select element,
    create a new DIV that will act as an option item:*/
    c = document.createElement("DIV");
    c.innerHTML = selElmnt.options[j].innerHTML;
    c.addEventListener("click", function(e) {
        /*when an item is clicked, update the original select box,
        and the selected item:*/
        var y, i, k, s, h, sl, yl;
        s = this.parentNode.parentNode.getElementsByTagName("select")[0];
        sl = s.length;
        h = this.parentNode.previousSibling;
        for (i = 0; i < sl; i++) {
          if (s.options[i].innerHTML == this.innerHTML) {
            s.selectedIndex = i;
            h.innerHTML = this.innerHTML;
            y = this.parentNode.getElementsByClassName("same-as-selected");
            yl = y.length;
            for (k = 0; k < yl; k++) {
              y[k].removeAttribute("class");
            }
            this.setAttribute("class", "same-as-selected");
            break;
          }
        }
        h.click();
    });
    b.appendChild(c);
  }
  x[i].appendChild(b);
  a.addEventListener("click", function(e) {
      /*when the select box is clicked, close any other select boxes,
      and open/close the current select box:*/
      e.stopPropagation();
      closeAllSelect(this);
      this.nextSibling.classList.toggle("select-hide");
      this.classList.toggle("select-arrow-active");
    });
}
function closeAllSelect(elmnt) {
  /*a function that will close all select boxes in the document,
  except the current select box:*/
  var x, y, i, xl, yl, arrNo = [];
  x = document.getElementsByClassName("select-items");
  y = document.getElementsByClassName("select-selected");
  xl = x.length;
  yl = y.length;
  for (i = 0; i < yl; i++) {
    if (elmnt == y[i]) {
      arrNo.push(i)
    } else {
      y[i].classList.remove("select-arrow-active");
    }
  }
  for (i = 0; i < xl; i++) {
    if (arrNo.indexOf(i)) {
      x[i].classList.add("select-hide");
    }
  }
}
/*if the user clicks anywhere outside the select box,
then close all select boxes:*/
document.addEventListener("click", closeAllSelect);
</script>

