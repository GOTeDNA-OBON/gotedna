window.onscroll = function() { scrollFunction(); };

function scrollFunction() {
  var scrollBtn = document.getElementById("scroll-top");
  if (!scrollBtn) return;  // exit if button doesn't exist

  if (document.documentElement.scrollTop > 300 || document.body.scrollTop > 300) {
    scrollBtn.style.display = "block";
  } else {
    scrollBtn.style.display = "none";
  }
}

function topFunction() {
  document.documentElement.scrollTop = 0;
  document.body.scrollTop = 0;
}

