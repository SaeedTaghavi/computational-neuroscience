// document.addEventListener("DOMContentLoaded", function () {
//   // search box placeholder
//   var input = document.querySelector(".md-search__input");
//   if (input) input.setAttribute("placeholder", "جستجو");

//   // "Type to start searching" — re-apply, since it renders after opening search
//   function fixSearchMeta() {
//     document.querySelectorAll(".md-search-result__meta").forEach(function (el) {
//       if (el.textContent.trim() === "Type to start searching") {
//         el.textContent = "برای جستجو تایپ کنید";
//       }
//     });
//   }
//   if (input) input.addEventListener("focus", function () {
//     setTimeout(fixSearchMeta, 50);
//   });
// });