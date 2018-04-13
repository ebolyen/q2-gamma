import loadWASM from './test.rs';

loadWASM().then(module => {
  const add = module.instance.exports['add'];
  console.log('foo:', add(1, 2));
})


window.LOAD_PACKED_TABLE = (payload) => {
  console.log("hello table")
}

window.LOAD_TREE = (payload) => {
  console.log("hello tree")
}
